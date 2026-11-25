import { useEffect, useRef, useState } from "react";
import { api } from "./api";
import type { FlowsheetModel, PropertyData, PropertyDataCheck } from "./types";

function checkBody(model: FlowsheetModel, data = model.property_data) {
  return JSON.stringify({
    components: model.components,
    property_method: model.property_method,
    custom_compounds: model.custom_compounds,
    henry_components: model.henry_components,
    property_data: data,
    units: model.units.filter((u) => u.property_method).map(({ id, type, property_method }) => ({ id, type, property_method })),
  });
}

export function PropertyDataEditor({ model, onChange }: {
  model: FlowsheetModel;
  onChange: (data: PropertyData) => void;
}) {
  const data = model.property_data;
  const body = checkBody(model);
  const hasComponents = model.components.length > 0;
  const text = JSON.stringify({ binary_parameters: data?.binary_parameters ?? [], unifac_groups: data?.unifac_groups ?? [] }, null, 2);
  const [draft, setDraft] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [saving, setSaving] = useState(false);
  const pendingApply = useRef<AbortController | null>(null);
  const [check, setCheck] = useState<{ body: string; report?: PropertyDataCheck; error?: string }>();
  const current = check?.body === body ? check : undefined;

  useEffect(() => {
    setDraft(null);
    setError(null);
  }, [text]);

  useEffect(() => () => pendingApply.current?.abort(), [body]);

  useEffect(() => {
    if (!hasComponents) return;
    const controller = new AbortController();
    const timer = window.setTimeout(async () => {
      try {
        const report = await api.checkPropertyData(body, controller.signal);
        if (!controller.signal.aborted) setCheck({ body, report });
      } catch (e) {
        if (!controller.signal.aborted) setCheck({ body, error: (e as Error).message });
      }
    }, 300);
    return () => {
      window.clearTimeout(timer);
      controller.abort();
    };
  }, [body, hasComponents]);

  async function apply() {
    const controller = new AbortController();
    pendingApply.current = controller;
    setError(null);
    setSaving(true);
    try {
      const parsed: unknown = JSON.parse(draft ?? text);
      if (!parsed || typeof parsed !== "object" || Array.isArray(parsed)) {
        throw new Error("Enter a JSON object with binary_parameters and unifac_groups.");
      }
      const next = { ...(parsed as PropertyData), strict: data?.strict ?? false };
      const nextBody = checkBody(model, next);
      const report = await api.checkPropertyData(nextBody, controller.signal);
      if (controller.signal.aborted) return;
      setCheck({ body: nextBody, report });
      onChange(next);
      setDraft(null);
    } catch (e) {
      if (!controller.signal.aborted) setError((e as Error).message);
    } finally {
      if (pendingApply.current === controller) {
        pendingApply.current = null;
        setSaving(false);
      }
    }
  }

  return (
    <>
      <label className="small">
        <input type="checkbox" checked={data?.strict ?? false} disabled={saving}
          onChange={(e) => onChange({ ...data, strict: e.target.checked })} /> Reject estimated or missing property data
      </label>
      <details>
        <summary>Local property data</summary>
        <p className="muted small">
          For each binary pair provide its method, components, convention (thermo), source, T_range_K,
          coefficients a and/or b in [ij, ji] order, and alpha for NRTL. Use P_range_Pa for a known pressure range.
          For groups, provide component, convention (dortmund), source, and groups as subgroup IDs with counts.
        </p>
        <textarea className="description" aria-label="Local property data JSON" spellCheck={false}
          rows={9} value={draft ?? text} disabled={saving} onChange={(e) => setDraft(e.target.value)} />
        <button onClick={apply} disabled={saving || draft === null}>{saving ? "Checking…" : "Apply data"}</button>
        {error && <p className="callout warn" role="alert">{error}</p>}
      </details>
      {hasComponents && (
        <div aria-live="polite">
          {!current && <p className="muted small">Checking property data…</p>}
          {current?.error && <p className="callout warn">Property data check failed. {current.error}</p>}
          {current?.report?.errors.map((err) => (
            <p className="callout warn" key={`${err.where}-${err.message}`}>
              {err.where && <b>{err.where} · </b>}{err.message} {err.hint}
            </p>
          ))}
          {current?.report?.sections.map((section) => (
            <div key={`${section.method}-${section.where}`}>
              <p className="small"><b>{section.where}</b> · {section.method}</p>
              {section.warnings.map((warning) => <p key={warning} className="callout warn">{warning}</p>)}
              {section.parameter_sources.length > 0 && (
                <details open={section.parameter_sources.some((s) => s.reference)}>
                  <summary>Parameter sources</summary>
                  <ul className="small">
                    {section.parameter_sources.map((s, i) => (
                      <li key={`${i}-${s.pair}`}>
                        <b>{s.pair}</b> · {s.source}{s.reference && <div>{s.reference}</div>}
                        {s.equation ? <div>{s.equation}</div> : s.convention && <div>{s.convention}</div>}
                        {s.T_range_K && <div>{s.T_range_K.join("–")} K</div>}
                        {s.P_range_Pa && <div>{s.P_range_Pa.join("–")} Pa</div>}
                      </li>
                    ))}
                  </ul>
                </details>
              )}
            </div>
          ))}
        </div>
      )}
    </>
  );
}
