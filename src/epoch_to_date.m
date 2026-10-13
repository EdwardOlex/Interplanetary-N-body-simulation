function s = epoch_to_date(epoch, t)
%EPOCH_TO_DATE  Calendar date of an elapsed time since the reference epoch.
%
%   s = EPOCH_TO_DATE(epoch, t) returns 'yyyy-mm-dd HH:MM' (TDB) for the
%   time T [s] elapsed since EPOCH (a datenum). Leap years are handled by
%   the calendar arithmetic of DATENUM. T may be a vector, in which case S
%   is a char matrix with one date per row. NaN gives '-'.

t = t(:);
s = repmat('-', numel(t), 16);
ok = ~isnan(t);
if any(ok)
    s(ok, :) = datestr(epoch + t(ok) / 86400, 'yyyy-mm-dd HH:MM');
end
s = strtrim(cellstr(s));
s = char(s);
end
