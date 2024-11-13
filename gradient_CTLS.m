% Fonction gradient de CTLS
function grad = gradient_CTLS(cx, cy)
    global xi yi R;
    n = length(xi);
    grad_cx = 0;
    grad_cy = 0;
    
    for i = 1:n
        Di = sqrt((xi(i) - cx)^2 + (yi(i) - cy)^2);
        if Di ~= 0
            grad_cx = grad_cx + 2 * (Di - R) * (cx - xi(i)) / Di;
            grad_cy = grad_cy + 2 * (Di - R) * (cy - yi(i)) / Di;
        end
    end
    
    grad = [grad_cx; grad_cy];
end