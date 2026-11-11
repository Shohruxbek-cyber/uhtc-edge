"""Download the NIST-JANAF tables this model reads and check them against the checksums of the tables used for the paper.

The tables are NIST Standard Reference Database 13 (doi 10.18434/T42S31). NIST reserves copyright on its Standard Reference
Data, so they are not redistributed with this repository. Each file is saved exactly as NIST serves it.

Usage: python literature/janaf/fetch_janaf.py
"""
import hashlib
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
URL = "https://janaf.nist.gov/tables/{}.txt"
TABLES = {  # file name: (NIST table id, SHA-256)
    "B2O3_cr": ("B-095", "e7576e70e9f8dfb408582a8e6c74b2203289859aaf0895da5f606edef5c4f30c"),
    "B2O3_cr-l": ("B-097", "bafed73b29d1161d1b5806f085c551d1183c8443b596544085fd815a6b7d5f11"),
    "B2O3_g": ("B-098", "80dc77f7444f7fd430af08e3b3f088846c88c95c6e3ca36a801c2054cee68705"),
    "B2O3_l": ("B-096", "d23fba24d5f85f6119d6836da13b4f43bdf3dc3d14e4311f4ac415f116be7295"),
    "CO2_g": ("C-095", "b9ecbb3758246a5d992e8313f281b81a7b33ddce20ae611dd5e8419e7b970aea"),
    "CO_g": ("C-093", "dcfea0f98f9ae41a91e7f4c1a713bec463a6f0cf52333216605a908ba73d9242"),
    "N2_ref": ("N-023", "e15d0f159ff4790511769075427bcdd3a3c04328209d40ddadd8726656f19621"),
    "NO_g": ("N-005", "f13273f065cfe9ab4b8e150750dcf65abbdb150834e9c763de63aef5b4a50587"),
    "N_g": ("N-002", "1619206bfaf1125f6531433856da5283c38385ed78d5c8f066721027db1012e5"),
    "O2_ref": ("O-029", "47bd5d4d600e962b0ca18c30016364762ce9d4ff1a07f39ff300648635f7e648"),
    "O_g": ("O-001", "cece32f5b417e5df12f22357934d4eeaa9cbacf751a8f76a17a62e41a57cca5e"),
    "SiC_cr_alpha": ("C-100", "1a079621dc763cd5e839ba00c8e80242327c0fc0e7ad2d086eb1c57d32de5d48"),
    "SiC_cr_beta": ("C-101", "91a489f7391c87aba1954b1e98189c5b1a2ea4ea90fbb978410c080e14d74627"),
    "SiO2_cr-l": ("O-039", "94d67d995f6dbe7eda10e4a0e43671cf6489e2a2fd58617bb064ca9f9df3512c"),
    "SiO2_cr_cristobalite-high": ("O-035", "b3b964ea6baf5ad47b27a212ef9192762c51fb2acf26f1cec521658fb41fb742"),
    "SiO2_cr_cristobalite-low": ("O-036", "6c98270832bd92260afeb266147ac7dc1fa073c355cddaa25d66803387ac2bd2"),
    "SiO2_cr_quartz": ("O-037", "9eabe66085680ba0bfb7a30c39a26c9ce5f45a4a463e6626b7ff0444b54b70b2"),
    "SiO2_g": ("O-040", "55f2f98221b348e89bdad573a36f83a0ed2bb7d4d3789b928a852e20a4a53845"),
    "SiO2_l": ("O-038", "f1233d320c6c49aa8055a95b7e0f2d6404b47d42079f42966c5f1331f94f3f2c"),
    "SiO_g": ("O-012", "e786b264b198a4e3fce25393c511dca8c24d2049651ff20c8f6f44079f064f45"),
    "ZrB2_cr": ("B-103", "d68dfa77f891ad827097093cb18f756ad9c941ca71a41b59a3b2cff538ba2cc2"),
    "ZrB2_cr-l": ("B-105", "2d1caed9768f56e61f62ce20da5482a324ae8e2255008db0574af77899b2da90"),
    "ZrB2_l": ("B-104", "77a09443089a8b197c37e4631e5ad1451a8f898fdb3bbf4e3499d535a1dad0e7"),
    "ZrO2_cr": ("O-052", "9db3b395bdcf5a6fdef36aed80ad653e4e0f1a5008cd39386292626047e1ba45"),
    "ZrO2_cr-l": ("O-054", "5eed109bbf57f72f069043b48e68f033ff6ae10499a5338e036817b682f11f35"),
    "ZrO2_l": ("O-053", "c13a6083aed93a2fe25c5d4d0352e9ce5bc1b5c8a8b6cfe59a97c99c9081c6ae"),
}


def main():
    bad = []
    for name, (nid, sha) in TABLES.items():
        f = HERE / f"{name}.txt"
        if f.exists() and hashlib.sha256(f.read_bytes()).hexdigest() == sha:
            continue
        req = urllib.request.Request(URL.format(nid), headers={"User-Agent": "uhtc-edge-setup/1.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read()
        if hashlib.sha256(data).hexdigest() != sha:
            bad.append(name)  # kept: NIST may have reformatted the file, and the user should know
        f.write_bytes(data)
        time.sleep(0.2)
    print(f"{len(TABLES)} JANAF tables in {HERE}; checksum differs from the published run for: {bad or 'none'}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
