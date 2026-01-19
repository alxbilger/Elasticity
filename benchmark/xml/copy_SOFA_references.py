import os
import shutil
import argparse


def main():
    parser = argparse.ArgumentParser(description='Copy regression tests from reference directory to target directory')
    parser.add_argument('ref_dir', help='Directory containing reference regression tests')
    parser.add_argument('target_dir', help='Directory to copy regression tests to')
    args = parser.parse_args()

    # Ensure target directory exists
    os.makedirs(args.target_dir, exist_ok=True)

    # Walk through reference directory
    for root, _, files in os.walk(args.ref_dir):
        for file in files:
            if file.endswith('.regression-tests') or file.endswith('json.gz'):
                # Calculate relative path from ref_dir
                rel_path = os.path.relpath(os.path.join(root, file), args.ref_dir)
                target_path = os.path.join(args.target_dir, rel_path)

                # Copy file while preserving path
                shutil.copy2(os.path.join(root, file), target_path)
                print(f"Copied: {os.path.join(root, file)} -> {target_path}")


if __name__ == "__main__":
    main()
