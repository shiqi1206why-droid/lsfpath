function field = apply_neumann_boundary(field)
%APPLY_NEUMANN_BOUNDARY Copy edge-adjacent values to ghost-cell boundaries.

    field(1, :) = field(2, :);
    field(end, :) = field(end-1, :);
    field(:, 1) = field(:, 2);
    field(:, end) = field(:, end-1);
end
