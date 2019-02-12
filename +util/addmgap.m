% ADDMGAP.
%
% MPF: I neglected to add the relative gap to the info_sol structure.
% This manually computes the gap for the natural image experiment. For
% `wflow`, just add 'NaN', since it's not relevant there and won't be
% used.

files = dir('cache/wusterland_saga*');
for i=1:length(files)
   name = files(i).name;
   file = sprintf('cache/%s',name);
   load(file);
   prObj = stats.info_sol.prObj;
   duObj = stats.info_sol.duObj;
   scale = stats.info_sol.scale;
   mGap = log10(prObj*duObj/scale);
   stats.info_sol.mGap = mGap;
   fprintf('%20s: %10.2e\n', name, mGap);
   save(file,'stats','solverOpts','genOpts');
end

files = dir('cache/wusterland_wflow*');
for i=1:length(files)
   name = files(i).name;
   file = sprintf('cache/%s',name);
   load(file);
   mGap = nan;
   stats.info_sol.mGap = mGap;
   fprintf('%20s: %10.2e\n', name, mGap);
   save(file,'stats','solverOpts','genOpts');
end

