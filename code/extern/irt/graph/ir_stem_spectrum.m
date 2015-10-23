function ir_stem_spectrum(x, y, thresh, varargin)
k = y > max(y(:)) * thresh;
stem(x(k), y(k), varargin{:})
