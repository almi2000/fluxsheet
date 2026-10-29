import type { Compound, Example, FlowsheetModel, Method, SolveResult, UnitType } from "./types";

async function request<T>(path: string, init?: RequestInit): Promise<T> {
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
  return res.json() as Promise<T>;
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
};
