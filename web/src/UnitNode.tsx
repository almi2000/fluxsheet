import { Handle, Position, type Node, type NodeProps } from "@xyflow/react";
import { memo } from "react";
import { getSymbol, portAt, UnitIcon } from "./icons";
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
  return (
    <div
      className={[
        "unit-node",
        selected ? "selected" : "",
        data.hasError ? "error" : "",
        data.stale ? "stale" : "",
      ].join(" ")}
      title={schema.description}
    >
      <div className="unit-symbol" style={{ width: sym.w, height: sym.h }}>
        <UnitIcon type={data.type} />
        {schema.inlets.map((p, i) => (
          <Handle
            key={p.name}
            id={p.name}
            type="target"
            position={Position.Left}
            style={at(portAt(sym, p.name, "in", i, schema.inlets.length))}
            title={p.label}
          >
            {schema.inlets.length > 1 && <span className="port-label in">{p.label}</span>}
          </Handle>
        ))}
        {schema.outlets.map((p, i) => (
          <Handle
            key={p.name}
            id={p.name}
            type="source"
            position={Position.Right}
            style={at(portAt(sym, p.name, "out", i, schema.outlets.length))}
            title={p.label}
          >
            {schema.outlets.length > 1 && <span className="port-label">{p.label}</span>}
          </Handle>
        ))}
      </div>
      <div className="unit-id">{id}</div>
      {badge && <div className="unit-badge">{badge}</div>}
      {data.method && <MethodTag method={data.method} override={!!data.property_method} />}
    </div>
  );
}

export const UnitNode = memo(UnitNodeView);
