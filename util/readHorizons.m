function data = readHorizons(fname,sv)
%INPUTS
%   fname   - (str) Full path to horizons output file
%   sv      - (bool) True if file contains velocity date. False if position
%              only (defaults false).
%    
%OUTPUTS
%   data     - Cell array containing each column from file
%
%NOTES
%
% Setup in Horizons is as follows:
%
% Ephemeris Type [change] : VECTORS
% Target Body [change] : 	Cassini (spacecraft) [-82]
% Coordinate Origin [change] : 	Saturn (body center) [500@699]
% Time Span [change] : 	Start=2006-01-01, Stop=2007-01-01, Step=1 d
% Table Settings [change] : quantities code=2; labels=NO; CSV format=YES
% Display/Output [change] : download/save (plain text file)

% Copyright (c) Dmitry Savransky 2020 (ds264@cornell.edu)

if ~exist(fname,'file')
    disp([fname,' does not exist.']);
    return
end
if ~exist('sv','var')
    sv = false;
end

fid = fopen(fname);
raw = textscan(fid, '%s', 'delimiter', '\n');
s = find(strcmpi(raw{1},'$$SOE'));
e = find(strcmpi(raw{1},'$$EOE'));
frewind(fid);
if sv
    fstring = '%f %s %f %f %f %f %f %f';
else
    fstring = '%f %s %f %f %f';
end
data = textscan(fid,fstring ,e-s+1,'delimiter', ',','headerlines',s);
fclose(fid);