"""Shared implementation for white-wine Tables S2 and S3.

Run either ``python code/TableS2_white_wine.py`` or
``python code/TableS3_white_wine.py`` from the repository root.

External Python dependencies: numpy, pandas, torch, scikit-learn.
The POM helper additionally requires R and the R package ``ordinal``.
"""

from __future__ import annotations

import copy
import math
import os
import pickle
import platform
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Sequence

import numpy as np
import pandas as pd
import torch
from sklearn import __version__ as sklearn_version
from sklearn.preprocessing import StandardScaler
from torch import nn
from torch.nn import functional as F


PROTOCOL_VERSION = "WHITE-WINE-TABLES-S2-S3-2026-08-28-v1"
MASTER_SEED = 1234
BATCH_SIZE = 32
MAX_EPOCHS = 500
PATIENCE = 30
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 1e-2
PROBABILITY_FLOOR = 1e-12

# The single 16-unit architecture is intentionally excluded.
ARCHITECTURES = (
    (16, 16),
    (16, 16, 16, 16),
    (32, 32),
    (64, 64),
)
BASELINES = ("S1", "S2", "S3", "S4", "S5")
PREDICTORS = (
    "fixed_acidity",
    "volatile_acidity",
    "citric_acid",
    "residual_sugar",
    "chlorides",
    "free_sulfur_dioxide",
    "total_sulfur_dioxide",
    "density",
    "ph",
    "sulphates",
    "alcohol",
)
MODEL_LABELS = {
    "cmnn": "CMNN GCLM",
    "softmax": "Cumulative softmax GCLM",
    "ci": "CI model",
    "si_cs": "SI-CS model",
    "pom": "POM",
}
TABLE_ORDER = tuple(MODEL_LABELS.values())


@dataclass(frozen=True)
class TableSpec:
    table_name: str
    n_classes: int
    response_map: dict[int, int]


@dataclass(frozen=True)
class Candidate:
    candidate_id: int
    architecture: tuple[int, ...]
    baseline: str | None = None


def configure_runtime() -> None:
    """Use the fixed CPU protocol used for the reported analysis."""

    torch.set_num_threads(1)
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass
    np.random.seed(MASTER_SEED)
    torch.manual_seed(MASTER_SEED)


def reset_seed() -> None:
    np.random.seed(MASTER_SEED)
    torch.manual_seed(MASTER_SEED)


def dependency_string() -> str:
    return (
        f"Python {platform.python_version()}; "
        f"NumPy {np.__version__}; pandas {pd.__version__}; "
        f"PyTorch {torch.__version__}; scikit-learn {sklearn_version}"
    )


def check_r_dependency(project_root: Path) -> None:
    """Fail before neural fitting if the POM software is unavailable."""

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise RuntimeError("Rscript is required but was not found on PATH")
    check = subprocess.run(
        [
            rscript,
            "-e",
            (
                "if (!requireNamespace('ordinal', quietly=TRUE)) "
                "quit(status=1) else cat(as.character(packageVersion('ordinal')))"
            ),
        ],
        cwd=project_root,
        capture_output=True,
        text=True,
    )
    if check.returncode != 0:
        raise RuntimeError(
            "R package 'ordinal' is required. Install it once with "
            "install.packages('ordinal')."
        )
    print(f"R ordinal {check.stdout.strip()}")


def _elu_mlp(
    n_inputs: int,
    hidden_widths: Sequence[int],
    n_outputs: int,
) -> nn.Sequential:
    layers: list[nn.Module] = []
    width_in = n_inputs
    for width_out in hidden_widths:
        layers.extend((nn.Linear(width_in, width_out), nn.ELU()))
        width_in = width_out
    layers.append(nn.Linear(width_in, n_outputs))
    return nn.Sequential(*layers)


class SIComplexShiftONTRAM(nn.Module):
    """Simple intercepts with one neural complex shift."""

    def __init__(
        self,
        n_predictors: int,
        n_classes: int,
        hidden_widths: Sequence[int],
    ) -> None:
        super().__init__()
        self.n_predictors = n_predictors
        self.n_classes = n_classes
        self.raw_gamma = nn.Parameter(torch.zeros(n_classes - 1))
        self.network = _elu_mlp(n_predictors, hidden_widths, 1)

    def transformation(self, x: torch.Tensor) -> torch.Tensor:
        theta = gamma_to_theta(self.raw_gamma)
        shift = self.network(x).squeeze(1)
        return theta.unsqueeze(0) - shift.unsqueeze(1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return log_probabilities_from_transformation(self.transformation(x))


class ComplexInterceptONTRAM(nn.Module):
    """Observation-specific ordered intercepts and zero shift."""

    def __init__(
        self,
        n_predictors: int,
        n_classes: int,
        hidden_widths: Sequence[int],
    ) -> None:
        super().__init__()
        self.n_predictors = n_predictors
        self.n_classes = n_classes
        self.network = _elu_mlp(n_predictors, hidden_widths, n_classes - 1)

    def transformation(self, x: torch.Tensor) -> torch.Tensor:
        return gamma_to_theta(self.network(x))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return log_probabilities_from_transformation(self.transformation(x))


def gamma_to_theta(raw_gamma: torch.Tensor) -> torch.Tensor:
    first = raw_gamma[..., :1]
    if raw_gamma.shape[-1] == 1:
        return first
    return torch.cat(
        (first, first + torch.cumsum(torch.exp(raw_gamma[..., 1:]), dim=-1)),
        dim=-1,
    )


def _log1mexp(nonpositive: torch.Tensor) -> torch.Tensor:
    split = -math.log(2.0)
    return torch.where(
        nonpositive < split,
        torch.log1p(-torch.exp(nonpositive)),
        torch.log(-torch.expm1(nonpositive)),
    )


def log_probabilities_from_transformation(h: torch.Tensor) -> torch.Tensor:
    if not torch.isfinite(h).all():
        raise FloatingPointError("The ONTRAM transformation is non-finite")
    if h.shape[1] > 1 and torch.any(h[:, 1:] - h[:, :-1] <= 0):
        raise FloatingPointError("The ONTRAM cutpoints are not strictly ordered")
    first = F.logsigmoid(h[:, :1])
    last = F.logsigmoid(-h[:, -1:])
    if h.shape[1] == 1:
        output = torch.cat((first, last), dim=1)
    else:
        lower = h[:, :-1]
        upper = h[:, 1:]
        middle = (
            F.logsigmoid(upper)
            + F.logsigmoid(-lower)
            + _log1mexp(lower - upper)
        )
        output = torch.cat((first, middle, last), dim=1)
    if not torch.isfinite(output).all():
        raise FloatingPointError("The ONTRAM log-probabilities are non-finite")
    return output


class PositiveIncrementRho(nn.Module):
    """Cumulative-softmax rho using one shared scalar ELU MLP."""

    def __init__(
        self,
        n_predictors: int,
        n_classes: int,
        hidden_widths: Sequence[int],
    ) -> None:
        super().__init__()
        self.n_predictors = n_predictors
        self.n_classes = n_classes
        layers: list[nn.Module] = []
        width_in = n_predictors + 1
        for width_out in hidden_widths:
            layers.extend((nn.Linear(width_in, width_out), nn.ELU()))
            width_in = width_out
        self.hidden = nn.Sequential(*layers)
        self.output = nn.Linear(width_in, 1)

    def components(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        n = x.shape[0]
        j = torch.arange(
            1, self.n_classes, dtype=x.dtype, device=x.device
        ) / self.n_classes
        x_grid = x[:, None, :].expand(n, self.n_classes - 1, self.n_predictors)
        j_grid = j[None, :, None].expand(n, self.n_classes - 1, 1)
        xj = torch.cat((x_grid, j_grid), dim=-1).reshape(
            n * (self.n_classes - 1), self.n_predictors + 1
        )
        relative_logits = self.output(self.hidden(xj)).reshape(
            n, self.n_classes - 1
        )
        logits = torch.cat(
            (relative_logits, torch.zeros((n, 1), device=x.device, dtype=x.dtype)),
            dim=1,
        )
        increments = torch.softmax(logits.to(torch.float64), dim=1)
        rho = self.n_classes * torch.cumsum(increments, dim=1)[:, :-1]
        return increments, rho


ACTIVATION_WEIGHTS = (7.0, 7.0, 2.0)


def combined_elu(x: torch.Tensor) -> torch.Tensor:
    units = x.shape[-1]
    total = sum(ACTIVATION_WEIGHTS)
    n_convex = round(ACTIVATION_WEIGHTS[0] / total * units)
    n_concave = round(ACTIVATION_WEIGHTS[1] / total * units)
    convex_x = x[..., :n_convex]
    concave_x = x[..., n_convex : n_convex + n_concave]
    saturated_x = x[..., n_convex + n_concave :]
    pieces: list[torch.Tensor] = []
    if convex_x.shape[-1]:
        pieces.append(F.elu(convex_x))
    if concave_x.shape[-1]:
        pieces.append(-F.elu(-concave_x))
    if saturated_x.shape[-1]:
        elu_one = F.elu(torch.ones((), dtype=x.dtype, device=x.device))
        left = F.elu(saturated_x + 1.0) - elu_one
        right = -F.elu(-(saturated_x - 1.0)) + elu_one
        pieces.append(torch.where(saturated_x <= 0.0, left, right))
    return torch.cat(pieces, dim=-1)


class CMNNDense(nn.Module):
    def __init__(
        self,
        in_features: int,
        out_features: int,
        monotonicity: Sequence[int],
        use_activation: bool,
    ) -> None:
        super().__init__()
        self.weight = nn.Parameter(torch.empty(in_features, out_features))
        self.bias = nn.Parameter(torch.zeros(out_features))
        self.register_buffer(
            "indicator", torch.as_tensor(monotonicity, dtype=torch.float32)
        )
        self.use_activation = use_activation
        nn.init.xavier_uniform_(self.weight)

    def constrained_weight(self) -> torch.Tensor:
        indicator = self.indicator[:, None]
        weight = torch.where(indicator == 1, torch.abs(self.weight), self.weight)
        return torch.where(indicator == -1, -torch.abs(self.weight), weight)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        output = x @ self.constrained_weight() + self.bias
        return combined_elu(output) if self.use_activation else output


class CMNNRho(nn.Module):
    """Partially monotone network, constrained only in the j/m input."""

    def __init__(
        self,
        n_predictors: int,
        n_classes: int,
        hidden_widths: Sequence[int],
    ) -> None:
        super().__init__()
        self.n_predictors = n_predictors
        self.n_classes = n_classes
        layers: list[CMNNDense] = []
        width_in = n_predictors + 1
        for layer_number, width_out in enumerate(hidden_widths):
            monotonicity = (
                [0] * n_predictors + [1]
                if layer_number == 0
                else [1] * width_in
            )
            layers.append(
                CMNNDense(width_in, width_out, monotonicity, use_activation=True)
            )
            width_in = width_out
        self.hidden = nn.ModuleList(layers)
        self.output = CMNNDense(width_in, 1, [1] * width_in, use_activation=False)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        n = x.shape[0]
        j = torch.arange(
            1, self.n_classes, dtype=x.dtype, device=x.device
        ) / self.n_classes
        x_grid = x[:, None, :].expand(n, self.n_classes - 1, self.n_predictors)
        j_grid = j[None, :, None].expand(n, self.n_classes - 1, 1)
        hidden = torch.cat((x_grid, j_grid), dim=-1).reshape(
            n * (self.n_classes - 1), self.n_predictors + 1
        )
        for layer in self.hidden:
            hidden = layer(hidden)
        f_value = torch.sigmoid(self.output(hidden))
        return self.n_classes * f_value.reshape(n, self.n_classes - 1)


def baseline_survival(
    rho: torch.Tensor,
    n_classes: int,
    baseline: str,
) -> torch.Tensor:
    eps = torch.finfo(rho.dtype).eps
    u = (rho / float(n_classes)).clamp(eps, 1.0 - eps)
    one_minus_u = 1.0 - u
    if baseline == "S1":
        return one_minus_u
    if baseline == "S2":
        return one_minus_u.square()
    if baseline == "S3":
        return torch.sqrt(one_minus_u)
    if baseline == "S4":
        return one_minus_u.square() / (one_minus_u.square() + u.square())
    return torch.sigmoid(-0.5 * (torch.log(u) - torch.log1p(-u)))


def category_probabilities_from_survival(survival: torch.Tensor) -> torch.Tensor:
    n = survival.shape[0]
    extended = torch.cat(
        (
            torch.ones((n, 1), dtype=survival.dtype, device=survival.device),
            survival,
            torch.zeros((n, 1), dtype=survival.dtype, device=survival.device),
        ),
        dim=1,
    )
    return extended[:, :-1] - extended[:, 1:]


def stable_softmax_probabilities(
    increments: torch.Tensor,
    baseline: str,
) -> torch.Tensor:
    """Evaluate S1--S5 category differences without cancellation."""

    if torch.any(increments <= 0) or not torch.isfinite(increments).all():
        raise FloatingPointError("A softmax increment is invalid")
    n = increments.shape[0]
    reverse_tail = torch.flip(
        torch.cumsum(torch.flip(increments, dims=(1,)), dim=1), dims=(1,)
    )
    ones = torch.ones((n, 1), dtype=increments.dtype, device=increments.device)
    zeros = torch.zeros_like(ones)
    tail_before = torch.cat((ones, reverse_tail[:, 1:]), dim=1).clamp(0.0, 1.0)
    tail_after = torch.cat((tail_before[:, 1:], zeros), dim=1)
    left = torch.cumsum(increments, dim=1)
    head_before = torch.cat((zeros, left[:, :-1]), dim=1).clamp(0.0, 1.0)
    head_after = torch.cat((left[:, :-1], ones), dim=1).clamp(0.0, 1.0)
    if baseline == "S1":
        output = increments
    elif baseline == "S2":
        output = increments * (tail_before + tail_after)
    elif baseline == "S3":
        output = increments / (torch.sqrt(tail_before) + torch.sqrt(tail_after))
    elif baseline == "S4":
        before = tail_before.square() + head_before.square()
        after = tail_after.square() + head_after.square()
        factor = tail_before * head_after + tail_after * head_before
        output = increments * factor / (before * after)
    else:
        first = torch.sqrt(increments[:, :1]) / (
            torch.sqrt(increments[:, :1])
            + torch.sqrt(increments[:, 1:].sum(dim=1, keepdim=True))
        )
        last = torch.sqrt(increments[:, -1:]) / (
            torch.sqrt(increments[:, -1:])
            + torch.sqrt(increments[:, :-1].sum(dim=1, keepdim=True))
        )
        if increments.shape[1] == 2:
            output = torch.cat((first, last), dim=1)
        else:
            tb = tail_before[:, 1:-1]
            ta = tail_after[:, 1:-1]
            hb = head_before[:, 1:-1]
            ha = head_after[:, 1:-1]
            cross = torch.sqrt(tb * ha) + torch.sqrt(ta * hb)
            interior = increments[:, 1:-1] / (
                cross * (torch.sqrt(tb) + torch.sqrt(hb))
                * (torch.sqrt(ta) + torch.sqrt(ha))
            )
            output = torch.cat((first, interior, last), dim=1)
    if not torch.isfinite(output).all():
        raise FloatingPointError("A category probability is non-finite")
    return output


def build_model(
    model_type: str,
    n_predictors: int,
    n_classes: int,
    architecture: tuple[int, ...],
) -> nn.Module:
    if model_type == "ci":
        return ComplexInterceptONTRAM(n_predictors, n_classes, architecture)
    if model_type == "si_cs":
        return SIComplexShiftONTRAM(n_predictors, n_classes, architecture)
    if model_type == "softmax":
        return PositiveIncrementRho(n_predictors, n_classes, architecture)
    if model_type == "cmnn":
        return CMNNRho(n_predictors, n_classes, architecture)
    raise ValueError(f"Unknown model type: {model_type}")


def model_log_probabilities(
    model_type: str,
    model: nn.Module,
    x: torch.Tensor,
    n_classes: int,
    baseline: str | None,
) -> tuple[torch.Tensor, torch.Tensor]:
    if model_type in {"ci", "si_cs"}:
        log_probabilities = model(x)
        probabilities = torch.exp(log_probabilities)
    elif model_type == "softmax":
        if baseline is None:
            raise ValueError("The cumulative-softmax model requires a baseline")
        increments, _ = model.components(x)
        probabilities = stable_softmax_probabilities(increments, baseline)
        log_probabilities = torch.log(probabilities.clamp_min(PROBABILITY_FLOOR))
    else:
        if baseline is None:
            raise ValueError("The CMNN model requires a baseline")
        rho = model(x)
        survival = baseline_survival(rho, n_classes, baseline)
        probabilities = category_probabilities_from_survival(survival)
        log_probabilities = torch.log(probabilities.clamp_min(PROBABILITY_FLOOR))
    if probabilities.shape != (x.shape[0], n_classes):
        raise AssertionError("The category-probability array has the wrong shape")
    if not torch.isfinite(log_probabilities).all():
        raise FloatingPointError("A category log-probability is non-finite")
    return log_probabilities, probabilities


def make_optimizer(model: nn.Module) -> torch.optim.AdamW:
    decay: list[nn.Parameter] = []
    no_decay: list[nn.Parameter] = []
    for name, parameter in model.named_parameters():
        (decay if name.endswith("weight") else no_decay).append(parameter)
    return torch.optim.AdamW(
        (
            {"params": decay, "weight_decay": WEIGHT_DECAY},
            {"params": no_decay, "weight_decay": 0.0},
        ),
        lr=LEARNING_RATE,
    )


def mean_nll(log_probabilities: torch.Tensor, labels: torch.Tensor) -> torch.Tensor:
    return -log_probabilities.gather(1, labels[:, None]).mean()


def candidate_grid(model_type: str) -> list[Candidate]:
    if model_type in {"ci", "si_cs"}:
        grid = [
            Candidate(candidate_id=i, architecture=architecture)
            for i, architecture in enumerate(ARCHITECTURES, start=1)
        ]
        if len(grid) != 4:
            raise AssertionError("A comparison model must have four candidates")
        return grid
    grid = [
        Candidate(candidate_id=i, architecture=architecture, baseline=baseline)
        for i, (baseline, architecture) in enumerate(
            (
                (baseline, architecture)
                for baseline in BASELINES
                for architecture in ARCHITECTURES
            ),
            start=1,
        )
    ]
    if len(grid) != 20:
        raise AssertionError("A GCLM must have twenty candidates")
    return grid


def _state_copy(model: nn.Module) -> dict[str, torch.Tensor]:
    return {
        name: value.detach().cpu().clone()
        for name, value in model.state_dict().items()
    }


def fit_candidate(
    model_type: str,
    candidate: Candidate,
    x_train: torch.Tensor,
    y_train: torch.Tensor,
    x_validation: torch.Tensor,
    y_validation: torch.Tensor,
    n_classes: int,
) -> dict[str, object]:
    reset_seed()
    model = build_model(
        model_type, x_train.shape[1], n_classes, candidate.architecture
    ).cpu()
    optimizer = make_optimizer(model)
    generator = torch.Generator(device="cpu")
    generator.manual_seed(MASTER_SEED)
    best_validation_nll = math.inf
    best_epoch = 0
    best_state: dict[str, torch.Tensor] | None = None
    epochs_without_improvement = 0
    start_time = perf_counter()

    for epoch in range(1, MAX_EPOCHS + 1):
        try:
            model.train()
            permutation = torch.randperm(x_train.shape[0], generator=generator)
            for start in range(0, x_train.shape[0], BATCH_SIZE):
                index = permutation[start : start + BATCH_SIZE]
                optimizer.zero_grad(set_to_none=True)
                log_probability, _ = model_log_probabilities(
                    model_type,
                    model,
                    x_train[index],
                    n_classes,
                    candidate.baseline,
                )
                loss = mean_nll(log_probability, y_train[index])
                if not torch.isfinite(loss):
                    raise FloatingPointError("The training NLL is non-finite")
                loss.backward()
                if any(
                    parameter.grad is None
                    or not torch.isfinite(parameter.grad).all()
                    for parameter in model.parameters()
                    if parameter.requires_grad
                ):
                    raise FloatingPointError("A training gradient is invalid")
                optimizer.step()

            model.eval()
            with torch.no_grad():
                validation_log_probability, _ = model_log_probabilities(
                    model_type,
                    model,
                    x_validation,
                    n_classes,
                    candidate.baseline,
                )
                validation_nll = float(
                    mean_nll(validation_log_probability, y_validation).item()
                )
        except (FloatingPointError, OverflowError, ValueError) as error:
            if best_state is None:
                raise RuntimeError(
                    f"Candidate {candidate.candidate_id} failed before a valid epoch"
                ) from error
            break

        if validation_nll < best_validation_nll:
            best_validation_nll = validation_nll
            best_epoch = epoch
            best_state = _state_copy(model)
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1
        if epochs_without_improvement >= PATIENCE:
            break

    if best_state is None or best_epoch < 1:
        raise RuntimeError(f"Candidate {candidate.candidate_id} has no valid epoch")
    return {
        "candidate": candidate,
        "best_epoch": best_epoch,
        "validation_nll": best_validation_nll,
        "elapsed_seconds": perf_counter() - start_time,
    }


def fit_selected_model(
    model_type: str,
    candidate: Candidate,
    selected_epochs: int,
    x_train: torch.Tensor,
    y_train: torch.Tensor,
    n_classes: int,
) -> tuple[nn.Module, int]:
    reset_seed()
    model = build_model(
        model_type, x_train.shape[1], n_classes, candidate.architecture
    ).cpu()
    optimizer = make_optimizer(model)
    generator = torch.Generator(device="cpu")
    generator.manual_seed(MASTER_SEED)
    last_valid_state: dict[str, torch.Tensor] | None = None
    completed_epochs = 0

    for epoch in range(1, selected_epochs + 1):
        try:
            model.train()
            permutation = torch.randperm(x_train.shape[0], generator=generator)
            for start in range(0, x_train.shape[0], BATCH_SIZE):
                index = permutation[start : start + BATCH_SIZE]
                optimizer.zero_grad(set_to_none=True)
                log_probability, _ = model_log_probabilities(
                    model_type,
                    model,
                    x_train[index],
                    n_classes,
                    candidate.baseline,
                )
                loss = mean_nll(log_probability, y_train[index])
                if not torch.isfinite(loss):
                    raise FloatingPointError("The refit NLL is non-finite")
                loss.backward()
                if any(
                    parameter.grad is None
                    or not torch.isfinite(parameter.grad).all()
                    for parameter in model.parameters()
                    if parameter.requires_grad
                ):
                    raise FloatingPointError("A refit gradient is invalid")
                optimizer.step()
            model.eval()
            with torch.no_grad():
                full_log_probability, _ = model_log_probabilities(
                    model_type,
                    model,
                    x_train,
                    n_classes,
                    candidate.baseline,
                )
                full_nll = mean_nll(full_log_probability, y_train)
                if not torch.isfinite(full_nll):
                    raise FloatingPointError("The full-training NLL is non-finite")
        except (FloatingPointError, OverflowError, ValueError) as error:
            if last_valid_state is None:
                raise RuntimeError("The selected model failed in its first epoch") from error
            model.load_state_dict(last_valid_state)
            break
        last_valid_state = _state_copy(model)
        completed_epochs = epoch

    if last_valid_state is None:
        raise RuntimeError("The selected model has no valid refit state")
    model.load_state_dict(last_valid_state)
    model.eval()
    return model, completed_epochs


def score_model(
    model_type: str,
    model: nn.Module,
    x_test: np.ndarray,
    y_test: np.ndarray,
    n_classes: int,
    baseline: str | None,
) -> tuple[np.ndarray, np.ndarray]:
    x_tensor = torch.as_tensor(x_test, dtype=torch.float32)
    y_tensor = torch.as_tensor(y_test, dtype=torch.long)
    model.eval()
    with torch.no_grad():
        if model_type == "cmnn":
            if baseline is None:
                raise ValueError("CMNN requires a baseline")
            rho = model(x_tensor).detach().cpu().to(torch.float64)
            probabilities = category_probabilities_from_survival(
                baseline_survival(rho, n_classes, baseline)
            )
            log_probability = torch.log(
                probabilities.clamp_min(PROBABILITY_FLOOR)
            )
        else:
            log_probability, probabilities = model_log_probabilities(
                model_type, model, x_tensor, n_classes, baseline
            )
        if float(torch.max(torch.abs(probabilities.sum(dim=1) - 1.0))) > 2e-6:
            raise AssertionError("Test probabilities do not sum to one")
        log_score = -log_probability.gather(1, y_tensor[:, None]).squeeze(1)
        predicted_cdf = torch.cumsum(probabilities, dim=1)[:, :-1]
        cutpoints = torch.arange(n_classes - 1, device=y_tensor.device)[None, :]
        observed_cdf = (y_tensor[:, None] <= cutpoints).to(probabilities.dtype)
        # This is the ordinary, unnormalized RPS: no division by m - 1.
        rps = torch.sum((predicted_cdf - observed_cdf).square(), dim=1)
    return log_score.cpu().numpy(), rps.cpu().numpy()


def load_inputs(project_root: Path, spec: TableSpec) -> tuple[pd.DataFrame, pd.DataFrame]:
    data_path = project_root / "data" / "wine_white_processed.csv"
    split_path = project_root / "data" / "white_repeated_stratified_5x5.csv"
    if not data_path.exists() or not split_path.exists():
        raise FileNotFoundError(
            "Run from the repository root and place both input CSV files directly "
            "under data/."
        )
    data = pd.read_csv(data_path)
    splits = pd.read_csv(split_path)
    missing = [name for name in ("quality", *PREDICTORS) if name not in data]
    if missing:
        raise ValueError(f"The data file is missing columns: {missing}")
    if "row_index" not in data:
        data.insert(0, "row_index", np.arange(len(data), dtype=np.int64))
    data["row_index"] = pd.to_numeric(data["row_index"], errors="raise").astype(int)
    if data["row_index"].duplicated().any():
        raise ValueError("row_index is duplicated")
    if not np.isfinite(data[list(PREDICTORS)].to_numpy(dtype=float)).all():
        raise ValueError("A predictor is missing or non-finite")
    quality = pd.to_numeric(data["quality"], errors="raise").astype(int)
    if not set(quality.unique()).issubset(spec.response_map):
        raise ValueError("quality contains a value outside the response map")
    data["y_index"] = quality.map(spec.response_map).astype(int) - 1
    if set(data["y_index"].unique()) != set(range(spec.n_classes)):
        raise ValueError("The recoded response does not contain every category")

    required = {"row_index", "repeat_id", "fold_id", "split"}
    if not required.issubset(splits):
        raise ValueError(f"The split file is missing: {sorted(required - set(splits))}")
    splits = splits[list(required)].copy()
    splits["row_index"] = pd.to_numeric(
        splits["row_index"], errors="raise"
    ).astype(int)
    splits["repeat_id"] = pd.to_numeric(
        splits["repeat_id"], errors="raise"
    ).astype(int)
    splits["fold_id"] = pd.to_numeric(splits["fold_id"], errors="raise").astype(int)
    if set(splits["repeat_id"]) != set(range(1, 6)):
        raise ValueError("repeat_id must be 1,...,5")
    if set(splits["fold_id"]) != set(range(1, 6)):
        raise ValueError("fold_id must be 1,...,5")
    if set(splits["split"].astype(str)) != {"inner_fit", "validation", "test"}:
        raise ValueError("split must contain inner_fit, validation, and test")
    for repeat_id in range(1, 6):
        tested = splits.loc[
            (splits["repeat_id"] == repeat_id) & (splits["split"] == "test"),
            "row_index",
        ]
        if len(tested) != len(data) or tested.nunique() != len(data):
            raise ValueError(f"Repeat {repeat_id} does not test every row exactly once")
    return data, splits


def fold_frames(
    data: pd.DataFrame,
    splits: pd.DataFrame,
    repeat_id: int,
    fold_id: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    assignment = splits.loc[
        (splits["repeat_id"] == repeat_id) & (splits["fold_id"] == fold_id),
        ["row_index", "split"],
    ]
    if len(assignment) != len(data) or assignment["row_index"].duplicated().any():
        raise ValueError("A repeat/fold assignment must contain every row once")
    frame = assignment.merge(data, on="row_index", validate="one_to_one")
    inner = frame.loc[frame["split"] == "inner_fit"].sort_values("row_index")
    validation = frame.loc[frame["split"] == "validation"].sort_values("row_index")
    test = frame.loc[frame["split"] == "test"].sort_values("row_index")
    if inner.empty or validation.empty or test.empty:
        raise ValueError("A fold component is empty")
    return inner, validation, test


def tune_and_score_fold(
    model_type: str,
    inner: pd.DataFrame,
    validation: pd.DataFrame,
    test: pd.DataFrame,
    n_classes: int,
) -> tuple[np.ndarray, np.ndarray, Candidate, int]:
    inner_scaler = StandardScaler()
    x_inner = inner_scaler.fit_transform(inner[list(PREDICTORS)]).astype(np.float32)
    x_validation = inner_scaler.transform(validation[list(PREDICTORS)]).astype(
        np.float32
    )
    y_inner = inner["y_index"].to_numpy(dtype=np.int64)
    y_validation = validation["y_index"].to_numpy(dtype=np.int64)
    x_inner_tensor = torch.from_numpy(x_inner)
    y_inner_tensor = torch.from_numpy(y_inner)
    x_validation_tensor = torch.from_numpy(x_validation)
    y_validation_tensor = torch.from_numpy(y_validation)

    results: list[dict[str, object]] = []
    candidates = candidate_grid(model_type)
    for index, candidate in enumerate(candidates, start=1):
        result = fit_candidate(
            model_type,
            candidate,
            x_inner_tensor,
            y_inner_tensor,
            x_validation_tensor,
            y_validation_tensor,
            n_classes,
        )
        results.append(result)
        print(
            f"      candidate {index:02d}/{len(candidates)}: "
            f"arch={'x'.join(map(str, candidate.architecture))}"
            + (f", baseline={candidate.baseline}" if candidate.baseline else "")
            + f", validation NLL={float(result['validation_nll']):.6f}, "
            f"epoch={int(result['best_epoch'])}",
            flush=True,
        )
    winner = min(
        results,
        key=lambda result: (
            float(result["validation_nll"]),
            int(result["candidate"].candidate_id),
        ),
    )
    selected: Candidate = winner["candidate"]
    selected_epochs = int(winner["best_epoch"])

    outer = pd.concat((inner, validation), ignore_index=True)
    outer_scaler = StandardScaler()
    x_outer = outer_scaler.fit_transform(outer[list(PREDICTORS)]).astype(np.float32)
    x_test = outer_scaler.transform(test[list(PREDICTORS)]).astype(np.float32)
    y_outer = outer["y_index"].to_numpy(dtype=np.int64)
    y_test = test["y_index"].to_numpy(dtype=np.int64)
    selected_model, completed_epochs = fit_selected_model(
        model_type,
        selected,
        selected_epochs,
        torch.from_numpy(x_outer),
        torch.from_numpy(y_outer),
        n_classes,
    )
    log_score, rps = score_model(
        model_type,
        selected_model,
        x_test,
        y_test,
        n_classes,
        selected.baseline,
    )
    return log_score, rps, selected, completed_epochs


def _save_progress(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        pickle.dump(payload, stream, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(temporary, path)


def _load_progress(path: Path, spec: TableSpec) -> dict[str, object]:
    if not path.exists():
        return {
            "protocol": PROTOCOL_VERSION,
            "table": spec.table_name,
            "n_classes": spec.n_classes,
            "folds": {},
        }
    with path.open("rb") as stream:
        payload = pickle.load(stream)
    expected = (PROTOCOL_VERSION, spec.table_name, spec.n_classes)
    observed = (payload.get("protocol"), payload.get("table"), payload.get("n_classes"))
    if observed != expected:
        raise RuntimeError(
            f"The progress file uses a different protocol: {path}. "
            "Delete it before restarting this table."
        )
    return payload


def run_neural_models(
    project_root: Path,
    spec: TableSpec,
    data: pd.DataFrame,
    splits: pd.DataFrame,
) -> pd.DataFrame:
    progress_path = project_root / "results" / f".{spec.table_name}_progress.pkl"
    progress = _load_progress(progress_path, spec)
    completed: dict[tuple[str, int, int], dict[str, object]] = progress["folds"]

    for model_type in ("ci", "si_cs", "softmax", "cmnn"):
        print(f"\n{MODEL_LABELS[model_type]}", flush=True)
        for repeat_id in range(1, 6):
            for fold_id in range(1, 6):
                key = (model_type, repeat_id, fold_id)
                if key in completed:
                    print(
                        f"  repeat {repeat_id}, fold {fold_id}: resumed",
                        flush=True,
                    )
                    continue
                print(f"  repeat {repeat_id}, fold {fold_id}", flush=True)
                inner, validation, test = fold_frames(
                    data, splits, repeat_id, fold_id
                )
                log_score, rps, selected, epochs = tune_and_score_fold(
                    model_type, inner, validation, test, spec.n_classes
                )
                completed[key] = {
                    "row_index": test["row_index"].to_numpy(dtype=np.int64),
                    "log_score": log_score,
                    "rps": rps,
                    "architecture": selected.architecture,
                    "baseline": selected.baseline,
                    "refit_epochs": epochs,
                }
                _save_progress(progress_path, progress)
                print(
                    f"    test LogS={log_score.mean():.6f}, "
                    f"RPS={rps.mean():.6f}",
                    flush=True,
                )

    repeat_rows: list[dict[str, object]] = []
    for model_type in ("ci", "si_cs", "softmax", "cmnn"):
        for repeat_id in range(1, 6):
            fold_parts = [
                completed[(model_type, repeat_id, fold_id)]
                for fold_id in range(1, 6)
            ]
            row_index = np.concatenate([part["row_index"] for part in fold_parts])
            if len(row_index) != len(data) or len(np.unique(row_index)) != len(data):
                raise AssertionError("A repeat does not contain one score per observation")
            repeat_rows.append(
                {
                    "model": MODEL_LABELS[model_type],
                    "repeat_id": repeat_id,
                    "LogS": float(
                        np.concatenate([part["log_score"] for part in fold_parts]).mean()
                    ),
                    "RPS": float(
                        np.concatenate([part["rps"] for part in fold_parts]).mean()
                    ),
                }
            )
    return pd.DataFrame(repeat_rows)


def run_pom(project_root: Path, spec: TableSpec) -> pd.DataFrame:
    helper = project_root / "code" / "white_wine_POM.R"
    if not helper.exists():
        raise FileNotFoundError(f"Missing POM helper: {helper}")
    with tempfile.TemporaryDirectory(prefix="white_wine_pom_") as directory:
        output_path = Path(directory) / "pom_repeat_scores.csv"
        command = [
            "Rscript",
            str(helper),
            str(spec.n_classes),
            str(output_path),
        ]
        try:
            subprocess.run(command, cwd=project_root, check=True)
        except FileNotFoundError as error:
            raise RuntimeError("Rscript is not installed or not on PATH") from error
        repeat_scores = pd.read_csv(output_path)
    if set(repeat_scores) != {"repeat_id", "LogS", "RPS"}:
        raise ValueError("The POM helper returned unexpected columns")
    if len(repeat_scores) != 5:
        raise ValueError("The POM helper must return five repeat scores")
    repeat_scores.insert(0, "model", MODEL_LABELS["pom"])
    return repeat_scores


def format_table(repeat_scores: pd.DataFrame) -> pd.DataFrame:
    summary = (
        repeat_scores.groupby("model", as_index=False)
        .agg(
            mean_LogS=("LogS", "mean"),
            sd_LogS=("LogS", "std"),
            mean_RPS=("RPS", "mean"),
            sd_RPS=("RPS", "std"),
        )
        .set_index("model")
        .reindex(TABLE_ORDER)
        .reset_index()
    )
    if summary.isna().any().any():
        raise AssertionError("The final table contains a missing value")
    return pd.DataFrame(
        {
            "Model": summary["model"],
            "LogS": [
                f"{mean:.4f} ({sd:.4f})"
                for mean, sd in zip(summary["mean_LogS"], summary["sd_LogS"])
            ],
            "RPS": [
                f"{mean:.4f} ({sd:.4f})"
                for mean, sd in zip(summary["mean_RPS"], summary["sd_RPS"])
            ],
        }
    )


def run_table(project_root: Path, spec: TableSpec) -> Path:
    project_root = project_root.resolve()
    configure_runtime()
    print(PROTOCOL_VERSION)
    print(dependency_string())
    print(f"Table: {spec.table_name}; response categories: {spec.n_classes}")
    print("Architectures: 2x16, 4x16, 2x32, 2x64")
    print("CI and SI-CS candidates per fold: 4")
    print("CMNN and cumulative-softmax candidates per fold: 20")
    print("RPS: ordinary unnormalized sum over cumulative cutpoints")

    data, splits = load_inputs(project_root, spec)
    check_r_dependency(project_root)
    neural_scores = run_neural_models(project_root, spec, data, splits)
    pom_scores = run_pom(project_root, spec)
    final_table = format_table(pd.concat((neural_scores, pom_scores), ignore_index=True))

    output_path = project_root / "results" / f"{spec.table_name}.csv"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.with_suffix(".csv.tmp")
    final_table.to_csv(temporary, index=False)
    os.replace(temporary, output_path)

    progress_path = project_root / "results" / f".{spec.table_name}_progress.pkl"
    progress_path.unlink(missing_ok=True)
    print(f"\n{spec.table_name}\n")
    print(final_table.to_string(index=False))
    print(f"\nSaved: {output_path}")
    return output_path
