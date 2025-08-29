function dist = cal_dist(p1, p2)
disp('dim-mm')
dist = sqrt(sum((p1-p2).^2))
end