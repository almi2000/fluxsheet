import { renderToStaticMarkup } from "react-dom/server";
import { describe, expect, it } from "vitest";
import { SolveStatus } from "./SolveStatus";
import type { SolveResult } from "./types";

const result: SolveResult = {
  ok: true, streams: {}, units: {}, property_method: null, property_method_reason: null,
  tear_streams: [], iterations: 1, converged: true, warnings: [], error: null,
};

describe("SolveStatus", () => {
  it("reports a missed specification separately from a closed balance", () => {
    const html = renderToStaticMarkup(<SolveStatus result={{
      ...result, specification_status: "warning", balance_status: "ok",
    }} />);
    expect(html).toContain("<b>Calculation</b> Converged");
    expect(html).toContain("<b>Duty specifications</b> Out of tolerance");
    expect(html).toContain("<b>Balances</b> Passed");
    expect(html).toContain("run-status warn");
  });

  it("shows unperformed checks after a calculation failure", () => {
    const html = renderToStaticMarkup(<SolveStatus result={{
      ...result, ok: false, converged: false, specification_status: "not_run", balance_status: "not_run",
    }} />);
    expect(html).toContain("<b>Calculation</b> Not completed");
    expect(html).toContain("<b>Duty specifications</b> Not run");
  });

  it("does not infer checks from older API responses", () => {
    expect(renderToStaticMarkup(<SolveStatus result={result} />)).toBe("");
  });
});
