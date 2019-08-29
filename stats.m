function Stats = stats(data)

Stats = struct('mean',mean(data),'std',std(data),...
    'min',min(data),'max',max(data),'median',median(data));