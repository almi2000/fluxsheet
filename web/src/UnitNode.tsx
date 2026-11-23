import { Handle, Position, type Node, type NodeProps } from "@xyflow/react";
import { memo } from "react";
import { getSymbol, portAt, UnitIcon } from "./icons";
import { activePorts } from "./mediumPorts";
import { MethodTag } from "./PropertyMethods";
import { formatQuantity } from "./quantities";
import type { Params, UnitResult, UnitType } from "./types";

export type UnitNodeData = {
  type: string;
  params: Params;
  schema: UnitType;
  result?: UnitResult;
  hasError?: boolean;
  stale?: boolean;
  property_method?: string; // the unit's own override, if any
  method?: string; // effective method, shown as a tag when the flowsheet mixes methods
  jacket?: "heat" | "cool"; // heated or cooled by a utility
  // Setup problems the editor sees before running (an idle side outlet) and the last run's
  // warnings about this unit, one per line.
  warning?: string;
};

export type UnitFlowNode = Node<UnitNodeData, "unit">;

/** The node badge: the unit's first energy result (duty or power), if non-zero. */
function summary(schema: UnitType, result?: UnitResult): string | null {
  if (!result) return null;
  const field = schema.result_fields.find(
    (f) => f.quantity === "power" && typeof result[f.name] === "number" && Math.abs(result[f.name] as number) > 1e-6,
  );
  return field ? formatQuantity(result[field.name], field.quantity) : null;
}

/** Handles sit on the symbol's nozzles; React Flow's own left/right placement is overridden. */
const at = (p: { x: number; y: number }) => ({ left: p.x, top: p.y, right: "auto", transform: "translate(-50%, -50%)" });

function UnitNodeView({ id, data, selected }: NodeProps<UnitFlowNode>) {
  const { schema } = data;
  const sym = getSymbol(data.type);
  const badge = summary(schema, data.result);
  const inlets = activePorts(schema, data.params, schema.inlets);
  const outlets = activePorts(schema, data.params, schema.outlets);
  return (
    <div
      className={[
        "unit-node",
        selected ? "selected" : "",
        data.hasError ? "error" : "",
        data.stale ? "stale" : "",
        data.warning ? "warning" : "",
      ].join(" ")}
      title={data.warning ? `${data.warning}

${schema.description}` : schema.description}
    >
      <div className="unit-symbol" style={{ width: sym.w, height: sym.h }}>
        <UnitIcon type={data.type} jacket={data.jacket} />
        {inlets.map((p, i) => (
          <Handle
            key={p.name}
            id={p.name}
            type="target"
            position={Position.Left}
            style={at(portAt(sym, p.name, "in", i, inlets.length))}
            title={p.label}
          >
            {inlets.length > 1 && <span className="port-label in">{p.label}</span>}
          </Handle>
        ))}
        {outlets.map((p, i) => (
          <Handle
            key={p.name}
            id={p.name}
            type="source"
            position={Position.Right}
            style={at(portAt(sym, p.name, "out", i, outlets.length))}
            title={p.label}
          >
            {outlets.length > 1 && <span className="port-label">{p.label}</span>}
          </Handle>
        ))}
      </div>
      <div className="unit-id">
        {id}
        {data.warning && (
          <span className="unit-warning" aria-label={data.warning}>
            {" "}⚠{data.warning.includes("\n") ? data.warning.split("\n").length : ""}
          </span>
        )}
      </div>
      {badge && <div className="unit-badge">{badge}</div>}
      {data.method && <MethodTag method={data.method} override={!!data.property_method} />}
    </div>
  );
}

export const UnitNode = memo(UnitNodeView);
