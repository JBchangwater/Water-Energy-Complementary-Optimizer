# WECO

**WECO** (Water-Energy Complementary Optimizer) is a physics-informed metaheuristic designed for general continuous optimization and complex engineering scheduling.

This repository provides a MATLAB implementation organized for direct GitHub release. The implementation follows the algorithmic logic reported in the accompanying paper, including:

- equivalent rainfall input constructed from inter-state differences,
- infiltration-runoff partitioning under water-energy complementarity,
- runoff-response routing with linear-reservoir dynamics,
- state-dependent regulation through soil moisture and hydraulic potential,
- feedback-based search using memory adaptation and a historical retention pool.

## Repository Structure

```text
WECO/
├── README.md
├── .gitignore
├── src/
│   └── WECO.m
├── examples/
│   ├── demo_sphere.m
│   └── demo_cec2017.m
└── benchmarks/
    └── cec2017/
        ├── Get_Functions_cec2017.m
        ├── cec17_func.mexw64
        └── input_data/
```

## Main Function

```matlab
[BestX, BestF, HisBestFit, out] = WECO(obj_fun, D, N0, MaxFEs, lb, ub)
[BestX, BestF, HisBestFit, out] = WECO(obj_fun, D, N0, MaxFEs, lb, ub, struct('seed', 2026))
```

### Inputs

- `obj_fun`: objective function handle.
- `D`: decision-space dimension.
- `N0`: initial population size.
- `MaxFEs`: maximum number of function evaluations.
- `lb`, `ub`: lower and upper bounds (scalar or 1-by-D vectors).
- `opts` (optional): either a struct or name-value pairs such as `seed`, `verbose`, and `history_length`.

### Outputs

- `BestX`: best solution found.
- `BestF`: best objective value found.
- `HisBestFit`: best-so-far history sampled in evaluation space.
- `out`: diagnostic structure containing run metadata and internal process statistics.

## Quick Start

### 1) Generic continuous optimization

Run the sphere example:

```matlab
run('examples/demo_sphere.m')
```

### 2) CEC2017 benchmark example

Run the benchmark example:

```matlab
run('examples/demo_cec2017.m')
```

## Notes

- The included `cec17_func.mexw64` file is a Windows MATLAB MEX binary.
- If you are working on another platform, replace it with a compatible CEC2017 implementation.
- The repository is intentionally lightweight and focuses on a clean reference implementation.

## Suggested Citation

If you use this code in your work, please cite the corresponding WECO paper.
