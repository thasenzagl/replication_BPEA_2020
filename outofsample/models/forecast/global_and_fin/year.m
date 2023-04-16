function out = year(time)
    time_vec = datevec(time);
    out = time_vec(:,1);
end 