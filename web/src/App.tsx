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
import { extractFragment, fragmentOrigin, parseFragment, remapFragment, type Point } from "./clipboard";
import { useHistory } from "./history";
import { FlowsheetInspector, StreamInspector, UnitInspector } from "./Inspector";
import { layeredLayout } from "./layout";
import { Palette } from "./Palette";
import { ResultsPanel } from "./ResultsPanel";
import { ShortcutHelp } from "./Shortcuts";
import { StreamEdge } from "./StreamEdge";
import type { Example, FlowsheetModel, Method, Params, Reaction, SolveResult, UnitType } from "./types";
import { UnitNode, type UnitFlowNode } from "./UnitNode";
import { fmt } from "./quantities";

const nodeTypes = { unit: UnitNode };
const edgeTypes = { stream: StreamEdge };

const ID_PREFIX: Record<string, string> = {
  Feed: "FEED", Product: "PROD", Mixer: "MIX", Splitter: "SPLIT", Heater: "HX", Flash: "FLASH",
  ComponentSplitter: "SEP", Valve: "VLV", Pump: "PUMP", Compressor: "COMP",
  StoichiometricReactor: "STOIC", YieldReactor: "YIELD", EquilibriumReactor: "EQ", GibbsReactor: "GIBBS",
  CSTR: "CSTR", PFR: "PFR",
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
};

type Meta = Pick<FlowsheetModel, "name" | "description" | "components" | "property_method"> & {
  reactions: Reaction[];
};
type Selection = { kind: "unit" | "stream"; id: string } | null;

const EMPTY_META: Meta = {
  name: "Untitled flowsheet", description: "", components: [], property_method: "auto", reactions: [],
};

function nextId(prefix: string, taken: Set<string>) {
  let i = 1;
  while (taken.has(`${prefix}${i}`)) i++;
  return `${prefix}${i}`;
}

function edgeStyle(e: Edge): Edge {
  return { ...e, type: "stream", markerEnd: { type: MarkerType.ArrowClosed, width: 16, height: 16 } };
}

/** Everything undo/redo restores. */
type Doc = { meta: Meta; nodes: UnitFlowNode[]; edges: Edge[] };

/** Canvas nodes and edges for a model's units and streams (used by Open and by Paste). */
function modelToCanvas(model: FlowsheetModel, schemas: Record<string, UnitType>) {
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
  const nodes: UnitFlowNode[] = units.map((u, i) => ({
    id: u.id,
    type: "unit",
    position: u.ui ?? { x: 60 + (i % 5) * 200, y: 80 + Math.floor(i / 5) * 180 },
    data: { type: u.type, params: u.params ?? {}, schema: schemas[u.type] },
  }));
  const edges = streams.map((s) => {
    const [source, sourceHandle] = s.from.split(".");
    const [target, targetHandle] = s.to!.split(".");
    return edgeStyle({ id: s.id, source, sourceHandle, target, targetHandle, label: s.id });
  });
  return { nodes, edges, positioned: units.every((u) => u.ui) };
}

/** Keyboard shortcuts stay out of the way while the user types in a field. */
function isTextField(t: EventTarget | null) {
  const el = t as HTMLElement | null;
  if (!el || !el.tagName) return false;
  if (el.isContentEditable || el.tagName === "TEXTAREA") return true;
  return el.tagName === "INPUT" && !["checkbox", "radio", "button", "range", "color", "file"].includes((el as HTMLInputElement).type);
}

const PASTE_OFFSET = 40;

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
  const [layoutOffer, setLayoutOffer] = useState(false);
  const [helpOpen, setHelpOpen] = useState(false);
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

  // ---- undo / redo ---------------------------------------------------------------------

  const docRef = useRef<Doc>({ meta, nodes, edges });
  docRef.current = { meta, nodes, edges };
  const currentDoc = useCallback(() => docRef.current, []);
  // Bumped on undo/redo so inspector fields with local drafts start over from the restored model.
  const [epoch, setEpoch] = useState(0);
  const applyDoc = useCallback(
    (d: Doc) => {
      setMeta(d.meta);
      setNodes(d.nodes);
      setEdges(d.edges);
      setSelection((s) =>
        s && (s.kind === "unit" ? d.nodes.some((n) => n.id === s.id) : d.edges.some((e) => e.id === s.id)) ? s : null,
      );
      setEpoch((n) => n + 1);
      setStale(true);
    },
    [setNodes, setEdges],
  );
  const history = useHistory(currentDoc, applyDoc);
  const [savedVersion, setSavedVersion] = useState(0);
  const dirty = history.version !== savedVersion;
  const dragStart = useRef<Doc | null>(null);

  useEffect(() => {
    if (!dirty) return;
    const warn = (e: BeforeUnloadEvent) => {
      e.preventDefault();
      e.returnValue = "";
    };
    window.addEventListener("beforeunload", warn);
    return () => window.removeEventListener("beforeunload", warn);
  }, [dirty]);

  // ---- model <-> canvas ----------------------------------------------------------------

  const toModel = useCallback(
    (): FlowsheetModel => ({
      format: "fluxsheet/1",
      ...meta,
      reactions: meta.reactions.length ? meta.reactions : undefined,
      units: nodes.map((n) => ({
        id: n.id,
        type: n.data.type,
        params: n.data.params,
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

  /**
   * Replaces the document (one undo step). `saved`: the new document counts as saved (Open);
   * otherwise it does only if there were no unsaved changes before (New, Examples).
   */
  const loadModel = useCallback(
    (model: FlowsheetModel, saved = false) => {
      const version = history.commit();
      if (saved || !dirty) setSavedVersion(version);
      const canvas = modelToCanvas(model, schemas);
      setNodes(canvas.nodes);
      setEdges(canvas.edges);
      setLayoutOffer(!canvas.positioned && canvas.nodes.length > 1);
      setMeta({
        name: model.name,
        description: model.description ?? "",
        components: model.components,
        property_method: model.property_method,
        reactions: model.reactions ?? [],
      });
      setResult(null);
      setStale(false);
      setSelection(null);
      setTimeout(() => flow.fitView({ padding: 0.2, duration: 300 }), 50);
    },
    [history, dirty, schemas, setNodes, setEdges, flow],
  );

  // Load the first example once the schemas are known, so a new user sees something work.
  const booted = useRef(false);
  useEffect(() => {
    if (booted.current || unitTypes.length === 0 || examples.length === 0) return;
    booted.current = true;
    api
      .example(examples[0].id)
      .then((m) => {
        loadModel(m, true);
        history.reset(); // the starting example is not an undo step
      })
      .catch(() => undefined);
  }, [unitTypes, examples, loadModel, history]);

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
      history.commit();
      setEdges((eds) => addEdge(edgeStyle({ ...c, id, label: id }), eds));
      touch();
    },
    [history, edges, setEdges, touch],
  );

  const onDrop = useCallback(
    (e: React.DragEvent) => {
      e.preventDefault();
      const type = e.dataTransfer.getData("application/fluxsheet-unit");
      if (!type || !schemas[type]) return;
      const position = flow.screenToFlowPosition({ x: e.clientX - 40, y: e.clientY - 30 });
      const id = nextId(ID_PREFIX[type] ?? type.toUpperCase(), new Set(nodes.map((n) => n.id)));
      history.commit();
      setNodes((ns) => [
        ...ns.map((n) => (n.selected ? { ...n, selected: false } : n)),
        { id, type: "unit", position, selected: true, data: { type, params: { ...(DEFAULT_PARAMS[type] ?? {}) }, schema: schemas[type] } },
      ]);
      setSelection({ kind: "unit", id });
      touch();
    },
    [history, flow, nodes, schemas, setNodes, touch],
  );

  const updateParams = (id: string, params: Params) => {
    history.commit(`params:${id}`);
    setNodes((ns) => ns.map((n) => (n.id === id ? { ...n, data: { ...n.data, params } } : n)));
    touch();
  };

  const renameUnit = (oldId: string, newId: string): string | null => {
    const clean = newId.trim();
    if (!clean || clean.includes(".")) return "Use letters, digits, - or _";
    if (nodes.some((n) => n.id === clean)) return `${clean} already exists`;
    history.commit();
    setNodes((ns) => ns.map((n) => (n.id === oldId ? { ...n, id: clean } : n)));
    setEdges((es) =>
      es.map((e) => ({ ...e, source: e.source === oldId ? clean : e.source, target: e.target === oldId ? clean : e.target })),
    );
    setSelection({ kind: "unit", id: clean });
    touch();
    return null;
  };

  const deleteUnits = (ids: ReadonlySet<string>, streamIds: ReadonlySet<string> = new Set()) => {
    if (ids.size === 0 && streamIds.size === 0) return;
    history.commit();
    setNodes((ns) => ns.filter((n) => !ids.has(n.id)));
    setEdges((es) => es.filter((e) => !streamIds.has(e.id) && !ids.has(e.source) && !ids.has(e.target)));
    setSelection(null);
    touch();
  };
  const deleteUnit = (id: string) => deleteUnits(new Set([id]));

  /** Units the next copy/delete acts on: the canvas selection, else the unit in the inspector. */
  const selectedUnitIds = (): Set<string> => {
    const ids = new Set(nodes.filter((n) => n.selected).map((n) => n.id));
    if (ids.size === 0 && selection?.kind === "unit") ids.add(selection.id);
    return ids;
  };

  const deleteSelection = () => {
    const streams = new Set(edges.filter((e) => e.selected).map((e) => e.id));
    if (selection?.kind === "stream") streams.add(selection.id);
    deleteUnits(selectedUnitIds(), streams);
  };

  // Dragging one or many units is one undo step, recorded when the drag ends with a move.
  const onDragStart = () => {
    dragStart.current ??= docRef.current;
  };
  const onDragStop = () => {
    const before = dragStart.current;
    dragStart.current = null;
    if (!before) return;
    const was = new Map(before.nodes.map((n) => [n.id, n.position]));
    const moved = docRef.current.nodes.some((n) => {
      const p = was.get(n.id);
      return p && (p.x !== n.position.x || p.y !== n.position.y);
    });
    if (moved) history.commitFrom(before);
  };

  // ---- copy / paste --------------------------------------------------------------------

  const clipboard = useRef<FlowsheetModel | null>(null);
  const pasteCount = useRef(0);
  const pointer = useRef<Point | null>(null); // screen position while over the canvas

  const copySelection = (): FlowsheetModel | null => {
    const ids = selectedUnitIds();
    if (ids.size === 0) return null;
    clipboard.current = extractFragment(toModel(), ids);
    pasteCount.current = 0;
    return clipboard.current;
  };

  /** Adds a fragment's units and inner streams with fresh ids at `at`; they become the selection. */
  const pasteFragment = (fragment: FlowsheetModel, at: Point) => {
    const known = { ...fragment, units: fragment.units.filter((u) => schemas[u.type]) };
    if (known.units.length === 0) return;
    const { units, streams } = remapFragment(known, new Set(nodes.map((n) => n.id)), new Set(edges.map((e) => e.id)), at);
    const canvas = modelToCanvas({ ...known, units, streams }, schemas);
    history.commit();
    setNodes((ns) => [
      ...ns.map((n) => (n.selected ? { ...n, selected: false } : n)),
      ...canvas.nodes.map((n) => ({ ...n, selected: true })),
    ]);
    setEdges((es) => [...es.map((e) => (e.selected ? { ...e, selected: false } : e)), ...canvas.edges]);
    // Units pasted from another flowsheet may refer to components this one lacks.
    const missing = (fragment.components ?? []).filter((c) => !meta.components.includes(c));
    if (missing.length) setMeta((m) => ({ ...m, components: [...m.components, ...missing] }));
    setSelection(units.length === 1 ? { kind: "unit", id: units[0].id } : null);
    touch();
  };

  /** Paste at the mouse pointer if it is over the canvas, else a bit below-right of the original. */
  const paste = (fragment: FlowsheetModel) => {
    if (pointer.current) return pasteFragment(fragment, flow.screenToFlowPosition(pointer.current));
    pasteCount.current += 1;
    const o = fragmentOrigin(fragment.units);
    const d = PASTE_OFFSET * pasteCount.current;
    pasteFragment(fragment, { x: o.x + d, y: o.y + d });
  };

  const duplicateSelection = () => {
    const ids = selectedUnitIds();
    if (ids.size === 0) return;
    const fragment = extractFragment(toModel(), ids);
    const o = fragmentOrigin(fragment.units);
    pasteFragment(fragment, { x: o.x + PASTE_OFFSET, y: o.y + PASTE_OFFSET });
  };

  const selectAll = () => {
    setNodes((ns) => ns.map((n) => ({ ...n, selected: true })));
    setSelection(null);
  };

  const clearSelection = () => {
    setNodes((ns) => ns.map((n) => (n.selected ? { ...n, selected: false } : n)));
    setEdges((es) => es.map((e) => (e.selected ? { ...e, selected: false } : e)));
    setSelection(null);
  };

  // ---- auto-layout ---------------------------------------------------------------------

  const arrange = () => {
    if (nodes.length === 0) return;
    const pos = layeredLayout(
      nodes.map((n) => ({
        id: n.id,
        width: n.measured?.width ?? 80,
        height: n.measured?.height ?? 70,
        role: n.data.type === "Feed" ? "feed" : n.data.type === "Product" ? "product" : "unit",
        inlets: n.data.schema?.inlets.map((p) => p.name) ?? [],
        outlets: n.data.schema?.outlets.map((p) => p.name) ?? [],
      })),
      edges,
    );
    history.commit();
    setNodes((ns) => ns.map((n) => ({ ...n, position: pos.get(n.id) ?? n.position })));
    setLayoutOffer(false);
    setTimeout(() => flow.fitView({ padding: 0.2, duration: 300 }), 50);
  };

  // ---- keyboard and clipboard events ---------------------------------------------------

  const onKeyDown = (e: KeyboardEvent) => {
    if (isTextField(e.target)) return;
    const mod = e.ctrlKey || e.metaKey;
    const key = e.key.toLowerCase();
    if (mod && key === "z") {
      if (e.shiftKey) history.redo();
      else history.undo();
    } else if (mod && key === "y") history.redo();
    else if (mod && (key === "c" || key === "x" || key === "v")) return clipboardKey(key);
    else if (mod && key === "d") duplicateSelection();
    else if (mod && key === "a") selectAll();
    else if (!mod && (e.key === "Delete" || e.key === "Backspace")) deleteSelection();
    else if (e.key === "?") setHelpOpen((o) => !o);
    else if (e.key === "Escape") {
      if (helpOpen) setHelpOpen(false);
      else clearSelection();
    } else return;
    e.preventDefault();
  };

  // Ctrl+C/X/V normally arrive as copy/cut/paste events, which can use the system clipboard.
  // Some browsers send no event when nothing editable has focus; then a fallback runs shortly
  // after the key press with the in-app clipboard (and the async clipboard API for copies).
  const pendingClipboard = useRef<number | undefined>(undefined);
  const clipboardKey = (key: string) => {
    // Leave ordinary text copying alone (selected text in the results table, ...).
    if (key !== "v" && window.getSelection()?.toString()) return;
    window.clearTimeout(pendingClipboard.current);
    if (key === "v") {
      pendingClipboard.current = window.setTimeout(() => clipboard.current && paste(clipboard.current), 50);
      return;
    }
    const fragment = copySelection();
    if (!fragment) return;
    pendingClipboard.current = window.setTimeout(() => {
      navigator.clipboard?.writeText(JSON.stringify(fragment, null, 2) + "\n").catch(() => undefined);
      if (key === "x") deleteSelection();
    }, 50);
  };

  const onCopy = (e: ClipboardEvent, cut: boolean) => {
    // Leave ordinary text copying alone (fields, selected text in the results table).
    if (isTextField(e.target) || window.getSelection()?.toString()) return;
    window.clearTimeout(pendingClipboard.current);
    const fragment = copySelection();
    if (!fragment) return;
    e.clipboardData?.setData("text/plain", JSON.stringify(fragment, null, 2) + "\n");
    e.preventDefault();
    if (cut) deleteSelection();
  };

  const onPaste = (e: ClipboardEvent) => {
    if (isTextField(e.target)) return;
    window.clearTimeout(pendingClipboard.current);
    const text = e.clipboardData?.getData("text/plain") ?? "";
    // The system clipboard wins (it may come from another tab); without text, use our own.
    const fragment = text ? parseFragment(text) : clipboard.current;
    if (!fragment) return;
    e.preventDefault();
    paste(fragment);
  };

  // The listeners are registered once and call the latest handlers through a ref.
  const handlers = useRef({ onKeyDown, onCopy, onPaste });
  handlers.current = { onKeyDown, onCopy, onPaste };
  useEffect(() => {
    const key = (e: KeyboardEvent) => handlers.current.onKeyDown(e);
    const copy = (e: ClipboardEvent) => handlers.current.onCopy(e, false);
    const cut = (e: ClipboardEvent) => handlers.current.onCopy(e, true);
    const paste = (e: ClipboardEvent) => handlers.current.onPaste(e);
    window.addEventListener("keydown", key);
    document.addEventListener("copy", copy);
    document.addEventListener("cut", cut);
    document.addEventListener("paste", paste);
    return () => {
      window.removeEventListener("keydown", key);
      document.removeEventListener("copy", copy);
      document.removeEventListener("cut", cut);
      document.removeEventListener("paste", paste);
    };
  }, []);

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
    setSavedVersion(history.version);
    history.seal();
  };

  const open = async (file: File) => {
    try {
      loadModel(JSON.parse(await file.text()), true);
    } catch (err) {
      alert(`Could not open ${file.name}: ${(err as Error).message}`);
    }
  };

  // ---- results onto the canvas ---------------------------------------------------------

  const errorAt = result && !result.ok ? result.error?.where ?? null : null;
  const shownNodes = useMemo(
    () =>
      nodes.map((n) => ({
        ...n,
        data: { ...n.data, result: result?.units[n.id], hasError: errorAt === n.id, stale },
      })),
    [nodes, result, errorAt, stale],
  );
  const shownEdges = useMemo(
    () =>
      edges.map((e) => {
        const s = result?.streams[e.id];
        return {
          ...e,
          selected: e.selected || (selection?.kind === "stream" && selection.id === e.id),
          className: [s && s.vapor_fraction >= 1 ? "vapor" : s && s.vapor_fraction > 0 ? "twophase" : "", errorAt === e.id ? "error" : ""].join(" "),
          label: s && !stale ? `${e.id} · ${fmt.T(s.T)} °C` : e.id,
        };
      }),
    [edges, result, selection, errorAt, stale],
  );

  const selectedNode = selection?.kind === "unit" ? nodes.find((n) => n.id === selection.id) : undefined;
  const selectedCount = nodes.reduce((k, n) => k + (n.selected ? 1 : 0), 0);
  const mod = /Mac|iPhone|iPad/.test(navigator.platform) ? "⌘" : "Ctrl+";

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
        <span className="doc-name" title={dirty ? "Unsaved changes" : undefined}>
          {meta.name}
          {dirty && <span className="dirty-dot" aria-label="unsaved changes" />}
        </span>
        <div className="toolbar-group" role="group" aria-label="Edit">
          <button className="icon-btn" onClick={history.undo} disabled={!history.canUndo} title={`Undo (${mod}Z)`} aria-label="Undo">
            <svg viewBox="0 0 16 16" aria-hidden="true"><path d="M5.5 3 2 6.5 5.5 10M2.5 6.5h7a4 4 0 0 1 0 8H7" /></svg>
          </button>
          <button className="icon-btn" onClick={history.redo} disabled={!history.canRedo} title={`Redo (${mod}Shift+Z)`} aria-label="Redo">
            <svg viewBox="0 0 16 16" aria-hidden="true"><path d="M10.5 3 14 6.5 10.5 10M13.5 6.5h-7a4 4 0 0 0 0 8H9" /></svg>
          </button>
          <button className="icon-btn" onClick={arrange} disabled={nodes.length === 0} title="Arrange the flowsheet left to right" aria-label="Arrange">
            <svg viewBox="0 0 16 16" aria-hidden="true">
              <rect x="1.5" y="6" width="3.5" height="4" rx="1" /><rect x="11" y="2" width="3.5" height="4" rx="1" />
              <rect x="11" y="10" width="3.5" height="4" rx="1" /><path d="M5 8h3m0-4v8m0-8h3m-3 8h3" />
            </svg>
            <span>Arrange</span>
          </button>
        </div>
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
        <button className="icon-btn" data-shortcut-toggle onClick={() => setHelpOpen((o) => !o)} title="Keyboard shortcuts (?)" aria-label="Keyboard shortcuts" aria-expanded={helpOpen}>
          ?
        </button>
        {helpOpen && <ShortcutHelp onClose={() => setHelpOpen(false)} />}
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
        {layoutOffer && (
          <div className="banner info">
            This flowsheet has no saved positions.
            <button className="primary" onClick={arrange}>Arrange automatically</button>
            <button onClick={() => setLayoutOffer(false)}>Keep as is</button>
          </div>
        )}
        {result?.ok && result.warnings.length > 0 && (
          <div className="banner warn">{result.warnings.map((w) => <div key={w}>⚠ {w}</div>)}</div>
        )}
      </div>

      <main className="workspace">
        <Palette types={unitTypes} />
        <div
          className="canvas"
          onDragOver={(e) => e.preventDefault()}
          onDrop={onDrop}
          onMouseMove={(e) => (pointer.current = { x: e.clientX, y: e.clientY })}
          onMouseLeave={() => (pointer.current = null)}
        >
          <ReactFlow
            nodes={shownNodes}
            edges={shownEdges}
            nodeTypes={nodeTypes}
            edgeTypes={edgeTypes}
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
            onNodeDragStart={onDragStart}
            onNodeDragStop={onDragStop}
            onSelectionDragStart={onDragStart}
            onSelectionDragStop={onDragStop}
            fitView
            snapToGrid
            snapGrid={[10, 10]}
            minZoom={0.2}
            deleteKeyCode={null}
            selectionKeyCode="Shift"
            multiSelectionKeyCode={["Shift", "Meta", "Control"]}
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
          {selectedCount > 1 ? (
            <div className="inspector-body">
              <h2>{selectedCount} units selected</h2>
              <p className="muted small">
                Drag to move them together, {mod}C / {mod}V to copy, {mod}D to duplicate.
              </p>
              <div className="multi-actions">
                <button onClick={duplicateSelection}>Duplicate</button>
                <button className="danger" onClick={deleteSelection}>Delete {selectedCount} units</button>
              </div>
            </div>
          ) : selectedNode ? (
            <UnitInspector
              key={epoch}
              id={selectedNode.id}
              schema={selectedNode.data.schema}
              params={selectedNode.data.params}
              components={meta.components}
              result={stale ? undefined : result?.units[selectedNode.id]}
              onRename={(newId) => renameUnit(selectedNode.id, newId)}
              onChange={(p) => updateParams(selectedNode.id, p)}
              onDelete={() => deleteUnit(selectedNode.id)}
            />
          ) : selection?.kind === "stream" ? (
            <StreamInspector id={selection.id} result={stale ? undefined : result?.streams[selection.id]} />
          ) : (
            <FlowsheetInspector
              key={epoch}
              name={meta.name}
              description={meta.description ?? ""}
              components={meta.components}
              method={meta.property_method}
              methods={methods}
              reactions={meta.reactions}
              onChange={(patch) => {
                // Typing a name or editing reactions coalesces into one step per field.
                history.commit(patch.components || patch.property_method ? null : `meta:${Object.keys(patch).join()}`);
                setMeta((m) => ({ ...m, ...patch }));
                if (patch.components || patch.property_method || patch.reactions) touch();
              }}
            />
          )}
        </aside>
      </main>

      <footer className="bottom">
        <ResultsPanel
          result={result}
          stale={stale}
          onSelectStream={(id) => setSelection({ kind: "stream", id })}
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
