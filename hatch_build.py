"""Bundle the built web editor and the example flowsheets into the wheel.

Both are optional for development installs (``pip install -e .``). For release builds set
``FLUXSHEET_REQUIRE_WEB=1`` so a wheel without the editor fails instead of shipping silently.
Build the editor first: ``cd web && npm ci && npm run build``.
"""

import os
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface

BUNDLED = {"web/dist": "fluxsheet/web", "examples": "fluxsheet/examples"}


class CustomBuildHook(BuildHookInterface):
    def initialize(self, version: str, build_data: dict) -> None:
        if version == "editable":
            return  # the source checkout already has both folders
        root = Path(self.root)
        for src, dest in BUNDLED.items():
            if (root / src).is_dir():
                build_data["force_include"][str(root / src)] = dest
            elif os.environ.get("FLUXSHEET_REQUIRE_WEB") == "1":
                raise FileNotFoundError(f"{src} is missing; build the web editor before releasing")
