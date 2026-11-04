import { useState } from "react";
import { api } from "./api";
import type { FlowsheetModel } from "./types";

/** Export buttons for the results: an Excel workbook, or a printable report in a new tab. */
export function ReportExport({ getModel }: { getModel: () => FlowsheetModel }) {
  const [busy, setBusy] = useState<"xlsx" | "html" | null>(null);

  const exportReport = async (format: "xlsx" | "html") => {
    // Open the tab inside the click, before the request, or popup blockers stop it.
    const tab = format === "html" ? window.open("", "_blank") : null;
    tab?.document.write("<p style='font-family:sans-serif'>Preparing the report…</p>");
    setBusy(format);
    try {
      const { blob, filename } = await api.report(getModel(), format);
      const url = URL.createObjectURL(blob);
      if (tab) {
        tab.location.href = url;
      } else {
        const a = document.createElement("a");
        a.href = url;
        a.download = filename;
        a.click();
      }
      setTimeout(() => URL.revokeObjectURL(url), 60_000); // the tab may still be loading it
    } catch (err) {
      tab?.close();
      alert(`Could not create the report: ${(err as Error).message}`);
    } finally {
      setBusy(null);
    }
  };

  return (
    <span className="report-export" title="Solves the current flowsheet and exports its results">
      Export:
      <button className="link" onClick={() => exportReport("xlsx")} disabled={busy !== null}>
        {busy === "xlsx" ? "Excel…" : "Excel"}
      </button>
      <button className="link" onClick={() => exportReport("html")} disabled={busy !== null}>
        {busy === "html" ? "Printable report…" : "Printable report"}
      </button>
    </span>
  );
}
