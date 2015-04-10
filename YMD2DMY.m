function DDMMMYYYY = YMD2DMY(YYYYMMDD)

if ~ischar(YYYYMMDD)
   error('Input date was of %s class, not char!', class(YYYYMMDD));
end
if numel(YYYYMMDD)~= 8
    error('Input date %s was not YYYYMMDD format!', YYYYMMDD);
end

YYYYMMDD = sprintf('%s/%s/%s', YYYYMMDD(1:4),YYYYMMDD(5:6),YYYYMMDD(7:8));
DDMMMYYYY = datestr(datenum(YYYYMMDD));