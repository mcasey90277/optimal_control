function v = optdef(s, f, dflt)
% OPTDEF  Optional struct-field default: return s.(f) if present AND nonempty,
% else dflt. Single source for the getdef/getf helpers formerly copy-pasted
% across the drivers.
% INPUTS:  s [struct]; f [char field name]; dflt [any default value]
% OUTPUTS: v - s.(f) if present and nonempty, else dflt
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
