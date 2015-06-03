function v = sr_clamp(v, vLowerB, vUpperB)

% SC_CLAMP: clamp value v

v = max(v, vLowerB);
v = min(v, vUpperB);

end