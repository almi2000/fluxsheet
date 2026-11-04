import type { Compound, DynamicResult, Example, FlowsheetModel, Method, OptimizationResult, SensitivityResult, SolveResult, UnitType } from "./types";

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const res = await send(path, init);
  return res.json() as Promise<T>;
}

/** fetch with JSON headers; a failed request throws with the server's explanation. */
async function send(path: string, init?: RequestInit): Promise<Response> {
  const res = await fetch(`/api${path}`, {
    headers: { "Content-Type": "application/json" },
    ...init,
  });
  if (!res.ok) {
    let detail = `${res.status} ${res.statusText}`;
    try {
      const body = await res.json();
      if (typeof body.detail === "string") detail = body.detail;
      else if (Array.isArray(body.detail))
        detail = body.detail.map((d: { msg: string }) => d.msg).join("; ");
    } catch {
      /* keep status text */
    }
    throw new Error(detail);
  }
  return res;
}

export const api = {
  unitTypes: () => request<UnitType[]>("/unit-types"),
  methods: () => request<Method[]>("/methods"),
  examples: () => request<Example[]>("/examples"),
  example: (id: string) => request<FlowsheetModel>(`/examples/${id}`),
  searchCompounds: (q: string) =>
    request<Compound[]>(`/compounds?q=${encodeURIComponent(q)}&limit=12`),
  suggestMethod: (components: string[]) =>
    request<{ method: string; reason: string }>("/methods/suggest", {
      method: "POST",
      body: JSON.stringify({ components }),
    }),
  solve: (model: FlowsheetModel) =>
    request<SolveResult>("/solve", { method: "POST", body: JSON.stringify(model) }),
  sensitivity: (model: FlowsheetModel, study: string) =>
    request<SensitivityResult>("/sensitivity", { method: "POST", body: JSON.stringify({ model, study }) }),
  optimize: (model: FlowsheetModel, optimization: string) =>
    request<OptimizationResult>("/optimize", { method: "POST", body: JSON.stringify({ model, optimization }) }),
  /** A report file of the solved flowsheet (see src/fluxsheet/report.py). */
  report: async (model: FlowsheetModel, format: "xlsx" | "html") => {
    const res = await send(`/report?format=${format}`, { method: "POST", body: JSON.stringify(model) });
    const name = /filename="([^"]+)"/.exec(res.headers.get("Content-Disposition") ?? "")?.[1];
    return { blob: await res.blob(), filename: name ?? `flowsheet_report.${format}` };
  },
  dynamic: (model: FlowsheetModel) =>
    request<DynamicResult>("/dynamic", { method: "POST", body: JSON.stringify(model) }),
};
