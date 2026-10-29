import { Handle, Position, type Node, type NodeProps } from "@xyflow/react";
import { memo } from "react";
import { UnitIcon } from "./icons";
import { formatQuantity } from "./quantities";
import type { Params, UnitResult, UnitType } from "./types";

export type UnitNodeData = {
  type: string;
  params: Params;
  schema: UnitType;
  result?: UnitResult;
  hasError?: boolean;
  stale?: boolean;
};

export type UnitFlowNode = Node<UnitNodeData, "unit">;

function spread(i: number, n: number) {
  return `${((i + 1) / (n + 1)) * 100}%`;
}

/** The node badge: the unit's first energy result (duty or power), if non-zero. */
function summary(schema: UnitType, result?: UnitResult): string | null {
  if (!result) return null;
  const field = schema.result_fields.find(
    (f) => f.quantity === "power" && typeof result[f.name] === "number" && Math.abs(result[f.name] as number) > 1e-6,
  );
  return field ? formatQuantity(result[field.name], field.quantity) : null;
}

function UnitNodeView({ id, data, selected }: NodeProps<UnitFlowNode>) {
  const { schema } = data;
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
      {schema.inlets.map((p, i) => (
        <Handle
          key={p.name}
          id={p.name}
          type="target"
          position={Position.Left}
          style={{ top: spread(i, schema.inlets.length) }}
          title={p.label}
        />
      ))}
      <UnitIcon type={data.type} />
      <div className="unit-id">{id}</div>
      {badge && <div className="unit-badge">{badge}</div>}
      {schema.outlets.map((p, i) => (
        <Handle
          key={p.name}
          id={p.name}
          type="source"
          position={Position.Right}
          style={{ top: spread(i, schema.outlets.length) }}
          title={p.label}
        >
          {schema.outlets.length > 1 && <span className="port-label">{p.label}</span>}
        </Handle>
      ))}
    </div>
  );
}

export const UnitNode = memo(UnitNodeView);
