#!/usr/bin/env python3
"""
Phase 0: Check local environment and dependencies.

Reports:
- Python version and platform
- CUDA availability
- Installed package versions
- System tools availability (mafft, cd-hit, mmseqs2)
- YAML config validity

Output: results/manifests/environment_report.json
"""

import json
import os
import sys
import platform
import subprocess
from pathlib import Path
from importlib.metadata import version as get_version

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


def get_python_info():
    """Return Python version and platform info."""
    return {
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "architecture": platform.machine(),
    }


def get_package_versions():
    """Return versions of key packages."""
    packages = [
        "biopython",
        "pandas",
        "numpy",
        "scipy",
        "scikit-learn",
        "umap-learn",
        "torch",
        "transformers",
        "py3Dmol",
        "matplotlib",
        "plotly",
        "seaborn",
        "tqdm",
        "pyyaml",
        "h5py",
    ]

    versions = {}
    for pkg in packages:
        try:
            # Map package names to import names
            import_name = {
                "umap-learn": "umap",
                "scikit-learn": "sklearn",
                "py3Dmol": "py3Dmol",
                "pyyaml": "yaml",
            }.get(pkg, pkg)
            versions[pkg] = get_version(pkg)
        except Exception as e:
            versions[pkg] = f"NOT INSTALLED ({e})"

    return versions


def check_cuda():
    """Check if CUDA is available."""
    try:
        import torch
        return {
            "cuda_available": torch.cuda.is_available(),
            "cuda_version": torch.version.cuda,
            "device_count": torch.cuda.device_count() if torch.cuda.is_available() else 0,
            "device_name": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        }
    except Exception as e:
        return {"cuda_error": str(e), "cuda_available": False}


def check_system_tools():
    """Check availability of optional system tools."""
    tools = ["mafft", "cd-hit", "mmseqs2"]
    results = {}

    for tool in tools:
        try:
            result = subprocess.run(
                [tool, "--version"],
                capture_output=True,
                timeout=5,
            )
            results[tool] = "available" if result.returncode == 0 else "not_found"
        except (FileNotFoundError, subprocess.TimeoutExpired):
            results[tool] = "not_found"

    return results


def check_config():
    """Validate config.yaml structure."""
    try:
        import yaml
        config_path = Path(__file__).parent.parent / "config" / "config.yaml"
        with open(config_path) as f:
            config = yaml.safe_load(f)

        required_sections = ["project", "ncbi", "taxa", "qc", "embeddings", "paths"]
        missing = [s for s in required_sections if s not in config]

        return {
            "config_found": True,
            "valid_yaml": True,
            "required_sections_present": len(missing) == 0,
            "missing_sections": missing,
            "config_path": str(config_path),
        }
    except FileNotFoundError:
        return {"config_found": False}
    except Exception as e:
        return {"config_found": True, "valid_yaml": False, "error": str(e)}


def main():
    report = {
        "timestamp": Path(
            __file__).parent.parent.joinpath("logs").mkdir(parents=True, exist_ok=True) or None,
        "python_info": get_python_info(),
        "packages": get_package_versions(),
        "cuda": check_cuda(),
        "system_tools": check_system_tools(),
        "config": check_config(),
        "summary": {},
    }

    # Generate summary
    missing_packages = [k for k, v in report["packages"].items() if "NOT INSTALLED" in v]
    report["summary"]["missing_packages"] = missing_packages if missing_packages else "all_installed"
    report["summary"]["cuda_available"] = report["cuda"].get("cuda_available", False)
    report["summary"]["required_system_tools_available"] = (
        report["system_tools"].get("mafft") == "available"
    )

    # Print to terminal
    print("\n" + "="*80)
    print("PHASE 0: ENVIRONMENT CHECK")
    print("="*80)
    print(f"\nPython: {report['python_info']['python_version']} ({report['python_info']['platform']})")
    print(f"CUDA available: {report['cuda'].get('cuda_available', False)}")
    if report['cuda'].get('cuda_available'):
        print(f"  Device: {report['cuda'].get('device_name')}")
        print(f"  CUDA version: {report['cuda'].get('cuda_version')}")

    print(f"\nMissing packages: {len(missing_packages)}")
    if missing_packages:
        for pkg in missing_packages:
            print(f"  - {pkg}")

    print(f"\nOptional tools:")
    for tool, status in report["system_tools"].items():
        print(f"  {tool}: {status}")

    print(f"\nConfig status: {report['config'].get('config_found', False)}")
    if report['config'].get('config_found'):
        print(f"  YAML valid: {report['config'].get('valid_yaml', False)}")
        print(f"  Required sections: {report['config'].get('required_sections_present', False)}")

    # Save report
    output_dir = Path(__file__).parent.parent / "results" / "manifests"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "environment_report.json"

    with open(output_file, "w") as f:
        json.dump(report, f, indent=2, default=str)

    print(f"\n✓ Report saved: {output_file}")
    print("="*80 + "\n")

    return 0 if not missing_packages else 1


if __name__ == "__main__":
    sys.exit(main())
