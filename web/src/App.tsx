import {
  Background,
  Controls,
  MarkerType,
  MiniMap,
  ReactFlow,
  ReactFlowProvider,
  addEdge,
  useEdgesState,
  useNodesState,
  useReactFlow,
  type Connection,
  type Edge,
  type IsValidConnection,
} from "@xyflow/react";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { api } from "./api";
import { FlowsheetInspector, StreamInspector, UnitInspector } from "./Inspector";
import { Palette } from "./Palette";
import { UnitMethodField } from "./PropertyMethods";
import { ResultsPanel } from "./ResultsPanel";
import { StudiesSettings } from "./StudiesEditor";
import type { DesignSpec, Example, FlowsheetModel, Method, Params, Reaction, SensitivityStudy, SolveResult, UnitType } from "./types";
import { UnitNode, type UnitFlowNode } from "./UnitNode";
import { fmt } from "./quantities";

const nodeTypes = { unit: UnitNode };

const ID_PREFIX: Record<string, string> = {
  Feed: "FEED", Product: "PROD", Mixer: "MIX", Splitter: "SPLIT", Heater: "HX", Flash: "FLASH",
  ComponentSplitter: "SEP", Valve: "VLV", Pump: "PUMP", Compressor: "COMP",
  StoichiometricReactor: "STOIC", YieldReactor: "YIELD", EquilibriumReactor: "EQ", GibbsReactor: "GIBBS",
  CSTR: "CSTR", PFR: "PFR", HeatExchanger: "EX",
};

const DEFAULT_PARAMS: Record<string, Params> = {
  Feed: { T: "25 degC", P: "1 bar", flows: {} },
  Heater: { T: "80 degC" },
  Valve: { P: "1 bar" },
  Pump: { P: "5 bar" },
  Compressor: { P: "5 bar" },
  ComponentSplitter: { fractions: {} },
  CSTR: { V: "1 m3", phase: "liquid" },
  PFR: { V: "1 m3", phase: "gas" },
  HeatExchanger: { min_approach: "10 K" },
};

type Meta = Pick<FlowsheetModel, "name" | "description" | "components" | "property_method"> & {
  reactions: Reaction[];
  design_specs: DesignSpec[];
  sensitivity: SensitivityStudy[];
};
type Selection = { kind: "unit" | "stream"; id: string } | null;

const EMPTY_META: Meta = {
  name: "Untitled flowsheet", description: "", components: [], property_method: "auto", reactions: [],
  design_specs: [], sensitivity: [],
};

function nextId(prefix: string, taken: Set<string>) {
  let i = 1;
  while (taken.has(`${prefix}${i}`)) i++;
  return `${prefix}${i}`;
}

function edgeStyle(e: Edge): Edge {
  return { ...e, type: "smoothstep", markerEnd: { type: MarkerType.ArrowClosed, width: 16, height: 16 } };
}

function Editor() {
  const [unitTypes, setUnitTypes] = useState<UnitType[]>([]);
  const [methods, setMethods] = useState<Method[]>([]);
  const [examples, setExamples] = useState<Example[]>([]);
  const [backendError, setBackendError] = useState<string | null>(null);

  const [meta, setMeta] = useState<Meta>(EMPTY_META);
  const [nodes, setNodes, onNodesChange] = useNodesState<UnitFlowNode>([]);
  const [edges, setEdges, onEdgesChange] = useEdgesState<Edge>([]);
  const [selection, setSelection] = useState<Selection>(null);
  const [result, setResult] = useState<SolveResult | null>(null);
  const [stale, setStale] = useState(false);
  const [running, setRunning] = useState(false);
  const fileInput = useRef<HTMLInputElement>(null);
  const flow = useReactFlow();

  const schemas = useMemo(() => Object.fromEntries(unitTypes.map((t) => [t.type, t])), [unitTypes]);

  useEffect(() => {
    Promise.all([api.unitTypes(), api.methods(), api.examples()])
      .then(([t, m, ex]) => {
        setUnitTypes(t);
        setMethods(m);
        setExamples(ex);
      })
      .catch(() =>
        setBackendError("Cannot reach the Fluxsheet server. Start it with: fluxsheet serve"),
      );
  }, []);

  const touch = useCallback(() => setStale(true), []);

  // ---- model <-> canvas ----------------------------------------------------------------

  const toModel = useCallback(
    (): FlowsheetModel => ({
      format: "fluxsheet/1",
      ...meta,
      reactions: meta.reactions.length ? meta.reactions : undefined,
      design_specs: meta.design_specs.length ? meta.design_specs : undefined,
      sensitivity: meta.sensitivity.length ? meta.sensitivity : undefined,
      units: nodes.map((n) => ({
        id: n.id,
        type: n.data.type,
        params: n.data.params,
        property_method: n.data.property_method,
        ui: { x: Math.round(n.position.x), y: Math.round(n.position.y) },
      })),
      streams: edges.map((e) => ({
        id: e.id,
        from: `${e.source}.${e.sourceHandle}`,
        to: `${e.target}.${e.targetHandle}`,
      })),
    }),
    [meta, nodes, edges],
  );

  const loadModel = useCallback(
    (model: FlowsheetModel) => {
      const units = [...model.units];
      const streams = model.streams.map((s) => ({ ...s }));
      // Streams leaving the flowsheet get an explicit Product so they show on the canvas.
      for (const s of streams) {
        if (!s.to) {
          const id = nextId("PROD", new Set(units.map((u) => u.id)));
          const src = units.find((u) => u.id === s.from.split(".")[0]);
          units.push({ id, type: "Product", ui: src?.ui ? { x: src.ui.x + 220, y: src.ui.y } : undefined });
          s.to = `${id}.in`;
        }
      }
      setNodes(
        units.map((u, i) => ({
          id: u.id,
          type: "unit",
          position: u.ui ?? { x: 60 + (i % 5) * 200, y: 80 + Math.floor(i / 5) * 180 },
          data: { type: u.type, params: u.params ?? {}, schema: schemas[u.type], property_method: u.property_method },
        })),
      );
      setEdges(
        streams.map((s) => {
          const [source, sourceHandle] = s.from.split(".");
          const [target, targetHandle] = s.to!.split(".");
          return edgeStyle({ id: s.id, source, sourceHandle, target, targetHandle, label: s.id });
        }),
      );
      setMeta({
        name: model.name,
        description: model.description ?? "",
        components: model.components,
        property_method: model.property_method,
        reactions: model.reactions ?? [],
        design_specs: model.design_specs ?? [],
        sensitivity: model.sensitivity ?? [],
      });
      setResult(null);
      setStale(false);
      setSelection(null);
      setTimeout(() => flow.fitView({ padding: 0.2, duration: 300 }), 50);
    },
    [schemas, setNodes, setEdges, flow],
  );

  // Load the first example once the schemas are known, so a new user sees something work.
  const booted = useRef(false);
  useEffect(() => {
    if (booted.current || unitTypes.length === 0 || examples.length === 0) return;
    booted.current = true;
    api.example(examples[0].id).then(loadModel).catch(() => undefined);
  }, [unitTypes, examples, loadModel]);

  // ---- editing -------------------------------------------------------------------------

  const isValidConnection: IsValidConnection = useCallback(
    (c) => {
      if (c.source === c.target) return false;
      const target = nodes.find((n) => n.id === c.target);
      const port = target?.data.schema.inlets.find((p) => p.name === c.targetHandle);
      const outletUsed = edges.some((e) => e.source === c.source && e.sourceHandle === c.sourceHandle);
      const inletUsed = edges.some((e) => e.target === c.target && e.targetHandle === c.targetHandle);
      return !outletUsed && (!inletUsed || !!port?.multiple);
    },
    [nodes, edges],
  );

  const onConnect = useCallback(
    (c: Connection) => {
      const id = nextId("S", new Set(edges.map((e) => e.id)));
      setEdges((eds) => addEdge(edgeStyle({ ...c, id, label: id }), eds));
      touch();
    },
    [edges, setEdges, touch],
  );

  const onDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      const type = e.dataTransfer.getData("application/fluxsheet-unit");
      if (!type || !schemas[type]) return;
      const position = flow.screenToFlowPosition({ x: e.clientX - 40, y: e.clientY - 30 });
      const id = nextId(ID_PREFIX[type] ?? type.toUpperCase(), new Set(nodes.map((n) => n.id)));
      setNodes((ns) => [
        ...ns,
        { id, type: "unit", position, data: { type, params: { ...(DEFAULT_PARAMS[type] ?? {}) }, schema: schemas[type] } },
      ]);
      setSelection({ kind: "unit", id });
      touch();
    },
    [flow, nodes, schemas, setNodes, touch],
  );

  const updateParams = (id: string, params: Params) => {
    setNodes((ns) => ns.map((n) => (n.id === id ? { ...n, data: { ...n.data, params } } : n)));
    touch();
  };

  const renameUnit = (oldId: string, newId: string): string | null => {
    const clean = newId.trim();
    if (!clean || clean.includes(".")) return "Use letters, digits, - or _";
    if (nodes.some((n) => n.id === clean)) return `${clean} already exists`;
    setNodes((ns) => ns.map((n) => (n.id === oldId ? { ...n, id: clean } : n)));
    setEdges((es) =>
      es.map((e) => ({ ...e, source: e.source === oldId ? clean : e.source, target: e.target === oldId ? clean : e.target })),
    );
    setSelection({ kind: "unit", id: clean });
    touch();
    return null;
  };

  const deleteUnit = (id: string) => {
    setNodes((ns) => ns.filter((n) => n.id !== id));
    setEdges((es) => es.filter((e) => e.source !== id && e.target !== id));
    setSelection(null);
    touch();
  };

  // ---- run / files ---------------------------------------------------------------------

  const run = async () => {
    setRunning(true);
    try {
      const r = await api.solve(toModel());
      setResult(r);
      setStale(false);
    } catch (err) {
      setResult({
        ok: false, streams: {}, units: {}, property_method: null, property_method_reason: null,
        tear_streams: [], iterations: 0, converged: false, warnings: [],
        error: { message: (err as Error).message, where: null, hint: null },
      });
    } finally {
      setRunning(false);
    }
  };

  const save = () => {
    const model = toModel();
    const blob = new Blob([JSON.stringify(model, null, 2) + "\n"], { type: "application/json" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${model.name.replace(/[^\w-]+/g, "_").toLowerCase() || "flowsheet"}.json`;
    a.click();
    URL.revokeObjectURL(a.href);
  };

  const open = async (file: File) => {
    try {
      loadModel(JSON.parse(await file.text()));
    } catch (err) {
      alert(`Could not open ${file.name}: ${(err as Error).message}`);
    }
  };

  // ---- results onto the canvas ---------------------------------------------------------

  const errorAt = result && !result.ok ? result.error?.where ?? null : null;
  // The flowsheet default as solved ("auto" resolved), and each unit's effective method.
  const defaultMethod = (!stale && result?.property_method) || meta.property_method;
  const unitMethod = useCallback(
    (n: UnitFlowNode): string | undefined =>
      (!stale ? (result?.units[n.id]?.property_method as string | undefined) : undefined) ??
      n.data.property_method ??
      (n.data.type === "Product" ? undefined : defaultMethod),
    [result, stale, defaultMethod],
  );
  const shownNodes = useMemo(() => {
    const mixed = new Set(nodes.map(unitMethod).filter(Boolean)).size > 1;
    return nodes.map((n) => ({
      ...n,
      data: {
        ...n.data, result: result?.units[n.id], hasError: errorAt === n.id, stale,
        method: mixed ? unitMethod(n) : undefined,
      },
    }));
  }, [nodes, result, errorAt, stale, unitMethod]);

  const setUnitMethod = (id: string, method: string | undefined) => {
    setNodes((ns) => ns.map((n) => (n.id === id ? { ...n, data: { ...n.data, property_method: method } } : n)));
    touch();
  };
  const shownEdges = useMemo(
    () =>
      edges.map((e) => {
        const s = result?.streams[e.id];
        return {
          ...e,
          selected: selection?.kind === "stream" && selection.id === e.id,
          className: [s && s.vapor_fraction >= 1 ? "vapor" : s && s.vapor_fraction > 0 ? "twophase" : "", errorAt === e.id ? "error" : ""].join(" "),
          label: s && !stale ? `${e.id} · ${fmt.T(s.T)} °C` : e.id,
        };
      }),
    [edges, result, selection, errorAt, stale],
  );

  const selectedNode = selection?.kind === "unit" ? nodes.find((n) => n.id === selection.id) : undefined;

  return (
    <div className="app">
      <header className="toolbar">
        <div className="brand">
          <svg width="22" height="22" viewBox="0 0 32 32" aria-hidden="true">
            <rect x="9" y="3" width="14" height="26" rx="7" fill="var(--accent)" />
            <path d="M9 16h14" stroke="var(--bg)" strokeWidth="2" />
          </svg>
          Fluxsheet
        </div>
        <span className="doc-name">{meta.name}</span>
        <div className="spacer" />
        <button onClick={() => loadModel({ format: "fluxsheet/1", ...EMPTY_META, units: [], streams: [] })}>New</button>
        <button onClick={() => fileInput.current?.click()}>Open…</button>
        <button onClick={save}>Save</button>
        <select
          className="examples"
          value=""
          onChange={(e) => e.target.value && api.example(e.target.value).then(loadModel)}
        >
          <option value="">Examples…</option>
          {examples.map((ex) => (
            <option key={ex.id} value={ex.id} title={ex.description}>
              {ex.name}
            </option>
          ))}
        </select>
        <button className="primary" onClick={run} disabled={running || !!backendError}>
          {running ? "Running…" : "▶ Run"}
        </button>
        <input
          ref={fileInput}
          type="file"
          accept=".json,application/json"
          hidden
          onChange={(e) => {
            const f = e.target.files?.[0];
            if (f) open(f);
            e.target.value = "";
          }}
        />
      </header>

      <div className="banners">
        {backendError && <div className="banner error">{backendError}</div>}
        {result && !result.ok && result.error && (
          <div className="banner error">
            <b>{result.error.where ? `${result.error.where}: ` : ""}</b>
            {result.error.message}
            {result.error.hint && <div className="hint">💡 {result.error.hint}</div>}
          </div>
        )}
        {result?.ok && result.warnings.length > 0 && (
          <div className="banner warn">{result.warnings.map((w) => <div key={w}>⚠ {w}</div>)}</div>
        )}
      </div>

      <main className="workspace">
        <Palette types={unitTypes} />
        <div className="canvas" onDragOver={(e) => e.preventDefault()} onDrop={onDrop}>
          <ReactFlow
            nodes={shownNodes}
            edges={shownEdges}
            nodeTypes={nodeTypes}
            onNodesChange={(c) => {
              onNodesChange(c);
              if (c.some((x) => x.type === "remove")) touch();
            }}
            onEdgesChange={(c) => {
              onEdgesChange(c);
              if (c.some((x) => x.type === "remove")) touch();
            }}
            onConnect={onConnect}
            isValidConnection={isValidConnection}
            onNodeClick={(_, n) => setSelection({ kind: "unit", id: n.id })}
            onEdgeClick={(_, e) => setSelection({ kind: "stream", id: e.id })}
            onPaneClick={() => setSelection(null)}
            fitView
            snapToGrid
            snapGrid={[10, 10]}
            deleteKeyCode={["Delete", "Backspace"]}
            proOptions={{ hideAttribution: true }}
          >
            <Background gap={20} size={1} />
            <Controls showInteractive={false} />
            <MiniMap pannable zoomable style={{ width: 150, height: 96 }} />
          </ReactFlow>
          {nodes.length === 0 && unitTypes.length > 0 && (
            <div className="canvas-hint">
              Drag a <b>Feed</b> from the left to start, or pick one of the <b>Examples</b>.
            </div>
          )}
        </div>
        <aside className="inspector">
          {selectedNode ? (
            <UnitInspector
              id={selectedNode.id}
              schema={selectedNode.data.schema}
              params={selectedNode.data.params}
              components={meta.components}
              result={stale ? undefined : result?.units[selectedNode.id]}
              onRename={(newId) => renameUnit(selectedNode.id, newId)}
              onChange={(p) => updateParams(selectedNode.id, p)}
              onDelete={() => deleteUnit(selectedNode.id)}
            >
              <UnitMethodField
                value={selectedNode.data.property_method}
                flowsheetMethod={defaultMethod}
                methods={methods}
                result={stale ? undefined : result?.units[selectedNode.id]}
                isProduct={selectedNode.data.type === "Product"}
                isSeparation={["Separators", "Columns"].includes(selectedNode.data.schema.category)}
                onChange={(m) => setUnitMethod(selectedNode.id, m)}
              />
            </UnitInspector>
          ) : selection?.kind === "stream" ? (
            <StreamInspector id={selection.id} result={stale ? undefined : result?.streams[selection.id]} />
          ) : (
            <FlowsheetInspector
              name={meta.name}
              description={meta.description ?? ""}
              components={meta.components}
              method={meta.property_method}
              methods={methods}
              reactions={meta.reactions}
              onChange={(patch) => {
                setMeta((m) => ({ ...m, ...patch }));
                if (patch.components || patch.property_method || patch.reactions) touch();
              }}
            >
              <StudiesSettings
                specs={meta.design_specs}
                studies={meta.sensitivity}
                results={stale ? undefined : result?.design_specs}
                units={nodes.map((n) => ({ id: n.id, schema: n.data.schema }))}
                streams={edges.map((e) => e.id)}
                components={meta.components}
                onChange={(patch) => {
                  setMeta((m) => ({ ...m, ...patch }));
                  if (patch.design_specs) touch();
                }}
              />
            </FlowsheetInspector>
          )}
        </aside>
      </main>

      <footer className="bottom">
        <ResultsPanel
          result={result}
          stale={stale}
          onSelectStream={(id) => setSelection({ kind: "stream", id })}
          studies={meta.sensitivity}
          getModel={toModel}
        />
      </footer>
    </div>
  );
}

export default function App() {
  return (
    <ReactFlowProvider>
      <Editor />
    </ReactFlowProvider>
  );
}
