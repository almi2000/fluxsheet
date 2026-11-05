import { useState } from "react";
import { api, DEMO_BODY_LIMIT_DETAIL } from "./api";
import { currentRunNames, currentRuns, useSettled, type Runs } from "./runs";
import type { FlowsheetModel } from "./types";

/**
 * Export buttons for the results: an Excel workbook, or a printable report in a new tab. The
 * report solves the flowsheet and adds the sensitivity studies, optimizations and (in Excel,
 * optionally) the dynamic run last run here, if the flowsheet has not changed since.
 */
export function ReportExport({
  getModel,
  getRuns,
  runCount = 0,
}: {
  getModel: () => FlowsheetModel;
  getRuns?: () => Runs;
  runCount?: number; // changes with every finished run
}) {
  const [busy, setBusy] = useState<"xlsx" | "html" | null>(null);
  const [withDynamics, setWithDynamics] = useState(true);
  // What would be included, for the tooltip; recomputed once edits settle, not per drag frame.
  const names = useSettled(() => (getRuns ? currentRunNames(getRuns(), getModel()) : []), [getModel, runCount]) ?? [];
  const hasDynamics = names.some((n) => n.startsWith("dynamic"));

  const exportReport = async (format: "xlsx" | "html") => {
    // Open the tab inside the click, before the request, or popup blockers stop it.
    const tab = format === "html" ? window.open("", "_blank") : null;
    tab?.document.write("<p style='font-family:sans-serif'>Preparing the report…</p>");
    setBusy(format);
    const m = getModel();
    // The time series only go to the workbook.
    const runs = getRuns ? currentRuns(getRuns(), m, withDynamics && format === "xlsx") : undefined;
    try {
      let file;
      try {
        file = await api.report(m, format, runs);
      } catch (err) {
        // The demo server takes small requests only: without the time series, else without runs.
        if (!runs || (err as Error).message !== DEMO_BODY_LIMIT_DETAIL) throw err;
        let left = "The runs are";
        try {
          if (!runs.dynamics) throw err;
          file = await api.report(m, format, { ...runs, dynamics: null });
          left = "The dynamic run is";
        } catch (again) {
          if ((again as Error).message !== DEMO_BODY_LIMIT_DETAIL) throw again;
          file = await api.report(m, format);
        }
        alert(`${left} too large for the demo server and left out of the report.`);
      }
      const url = URL.createObjectURL(file.blob);
      if (tab) {
        tab.location.href = url;
      } else {
        const a = document.createElement("a");
        a.href = url;
        a.download = file.filename;
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

  const title =
    "Solves the current flowsheet and exports its results" +
    (names.length ? `, with ${names.join(", ")}` : "") +
    ". Studies, optimizations and dynamic runs are included when run since the last change.";
  return (
    <span className="report-export" title={title}>
      Export:
      <button className="link" onClick={() => exportReport("xlsx")} disabled={busy !== null}>
        {busy === "xlsx" ? "Excel…" : "Excel"}
      </button>
      <button className="link" onClick={() => exportReport("html")} disabled={busy !== null}>
        {busy === "html" ? "Printable report…" : "Printable report"}
      </button>
      {hasDynamics && (
        <label className="small" title="Add the time series of the last dynamic run as a sheet of the Excel report">
          <input type="checkbox" checked={withDynamics} onChange={(e) => setWithDynamics(e.target.checked)} /> with
          dynamic run
        </label>
      )}
    </span>
  );
}
