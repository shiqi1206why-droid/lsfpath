function write_path_to_file(csv_path, xr, yr, xs, ys)
%WRITE_PATH_TO_FILE Export one raw/smoothed segment to CSV.

    n = max([numel(xr), numel(yr), numel(xs), numel(ys)]);
    raw_x = nan(n, 1);
    raw_y = nan(n, 1);
    smooth_x = nan(n, 1);
    smooth_y = nan(n, 1);
    raw_x(1:numel(xr)) = xr(:);
    raw_y(1:numel(yr)) = yr(:);
    smooth_x(1:numel(xs)) = xs(:);
    smooth_y(1:numel(ys)) = ys(:);
    T = table(raw_x, raw_y, smooth_x, smooth_y);
    writetable(T, csv_path);
end
