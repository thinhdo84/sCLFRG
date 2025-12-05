function us = sat_scalar(u,u_up,u_low)
us = min(u_up, max(u_low, u));
end