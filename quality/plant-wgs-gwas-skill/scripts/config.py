#!/usr/bin/env python3
"""读取 project_config.yaml 的轻量工具（仅依赖 PyYAML）。

用法：
  python3 config.py <config.yaml> get <dot.path>      # 取单个值，支持 phenotypes.0.name
  python3 config.py <config.yaml> count <dot.path>    # 取列表长度
  python3 config.py <config.yaml> validate            # 校验必填字段，输出 OK / 缺项
  python3 config.py <config.yaml> dump                # 打印完整配置（调试）
"""
import sys
import yaml

REQUIRED = [
    "project.name", "species.name", "species.reference_genome", "species.chr_count",
    "population.type", "paths.workdir", "genotype.vcf",
    "qc.mind", "qc.geno", "qc.maf", "ld_pruning.window", "ld_pruning.r2",
    "structure.n_pc", "phenotypes", "model.default", "threshold.method",
]


def load(path):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def resolve(cfg, dotted):
    cur = cfg
    for part in dotted.split("."):
        if isinstance(cur, list):
            cur = cur[int(part)]
        elif isinstance(cur, dict) and part in cur:
            cur = cur[part]
        else:
            raise KeyError(dotted)
    return cur


def scalar(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (list, dict)):
        return yaml.safe_dump(v, default_flow_style=True).strip()
    return str(v)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    cfg_path, action = sys.argv[1], sys.argv[2]
    key = sys.argv[3] if len(sys.argv) > 3 else ""
    cfg = load(cfg_path)
    if action == "get":
        print(scalar(resolve(cfg, key)))
    elif action == "count":
        v = resolve(cfg, key)
        print(len(v) if isinstance(v, list) else 0)
    elif action == "dump":
        print(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False))
    elif action == "validate":
        missing = []
        for k in REQUIRED:
            try:
                resolve(cfg, k)
            except (KeyError, IndexError, TypeError, ValueError):
                missing.append(k)
        if missing:
            print("MISSING: " + ", ".join(missing))
            sys.exit(1)
        print("OK")
    else:
        print("unknown action: " + action, file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
