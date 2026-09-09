"""Print a Plotly figure using existing Windows Chrome from WSL."""
from pathlib import Path
import subprocess
import tempfile
import shutil

def windows_path(path):
    return subprocess.check_output(["wslpath", "-w", str(path)], text=True).strip()

def export_pdf(fig, output, chrome):
    chrome = Path(chrome).resolve()
    if not chrome.is_file():
        raise ValueError(f"Chrome executable not found: {chrome}")
    # Chrome's usual Windows installation is under Local/Google/Chrome/Application.
    windows_temp = chrome.parents[3] / "Temp"
    if not windows_temp.is_dir():
        raise ValueError(f"Windows temporary directory not found: {windows_temp}")
    with tempfile.TemporaryDirectory(prefix="sequencing-plot-", dir=windows_temp) as tmp:
        tmp = Path(tmp)
        html = tmp / "figure.html"
        pdf = tmp / "figure.pdf"
        width, height = fig.layout.width, fig.layout.height
        content = fig.to_html(include_plotlyjs=True, full_html=True, config={"responsive": False, "displayModeBar": False})
        style = f"<style>@page {{size: {width}px {height}px; margin: 0;}} html,body {{margin:0;padding:0;width:{width}px;height:{height}px;overflow:hidden;}}</style>"
        html.write_text(content.replace("</head>", style + "</head>"), encoding="utf-8")
        uri = "file:///" + windows_path(html).replace("\\", "/")
        command = [str(chrome), "--headless=new", "--disable-gpu", "--no-first-run",
                   "--no-default-browser-check", "--no-pdf-header-footer",
                   "--user-data-dir=" + windows_path(tmp / "profile"),
                   "--print-to-pdf=" + windows_path(pdf), "--virtual-time-budget=10000",
                   "--run-all-compositor-stages-before-draw", uri]
        result = subprocess.run(command, capture_output=True, timeout=90)
        if result.returncode or not pdf.exists() or pdf.stat().st_size < 1000:
            raise RuntimeError("Windows Chrome PDF export failed: " + result.stderr.decode(errors="replace")[-3000:])
        shutil.copyfile(pdf, output)
    return {"backend": "Windows Chrome headless print of identical Plotly figure", "browser": str(chrome)}
