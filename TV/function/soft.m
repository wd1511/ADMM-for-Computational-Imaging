function u = soft(s, mu, gamma,sep)
sa = real(s); sb=imag(s);
if sep
    ua = sign(sa).*max(abs(sa)-mu/gamma,0);
    ub = sign(sb).*max(abs(sb)-mu/gamma,0);
else
    coef = sqrt(sa.^2./(sa.^2+sb.^2));
    ua = coef.*sign(sa).* max(sqrt(sa.^2+sb.^2)-mu/gamma,0);
    ub = sb ./ sa .* ua;
end
u = ua + ub * sqrt(-1);