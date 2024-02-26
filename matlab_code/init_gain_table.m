function GT = init_gain_table( varargin )
%INIT_GAIN_TABLE Creates an empty gain table
%   Optional parameters:
%   - MAX_CLIQUE_SIZE (default 6)
%   - INITIAL_GAIN_TABLE_SIZE (default 100)
if nargin == 2 
    GT.MAX_CLIQUE_SIZE = varargin{1};
    GT.INITIAL_GAIN_TABLE_SIZE = varargin{2};
elseif nargin == 1
    GT.MAX_CLIQUE_SIZE = varargin{1};
    GT.INITIAL_GAIN_TABLE_SIZE = 100;
else
    GT.MAX_CLIQUE_SIZE = 6;
    GT.INITIAL_GAIN_TABLE_SIZE = 100;
end;

GT.MAX_SEPARATOR_SIZE = GT.MAX_CLIQUE_SIZE - 1;
GT.rowid              = NaN(GT.INITIAL_GAIN_TABLE_SIZE, 1);
GT.cliques            = NaN(GT.INITIAL_GAIN_TABLE_SIZE, GT.MAX_CLIQUE_SIZE);
GT.separators         = NaN(GT.INITIAL_GAIN_TABLE_SIZE, GT.MAX_SEPARATOR_SIZE);
GT.coordination_num   = NaN(GT.INITIAL_GAIN_TABLE_SIZE, 1);
GT.gains              = NaN(GT.INITIAL_GAIN_TABLE_SIZE, 1);
GT.nodes              = NaN(GT.INITIAL_GAIN_TABLE_SIZE, 1);
GT.tot_records        = 0;
end

