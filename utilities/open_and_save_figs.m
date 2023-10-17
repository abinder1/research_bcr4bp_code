clear; clc; close('all')

listing = dir("..\figures\JGCD 2023\.fig\")
dest_folder = dir("..\figures\JGCD 2023\.eps\")

for k = 3:1:length(listing)
    filename = listing(k).name
    full_path = fullfile(listing(k).folder, filename)
    save_path = fullfile(dest_folder(1).folder, filename(1:end-4))

    f = openfig(full_path, 'invisible');
    saveas(f, save_path, 'epsc');
    close(f)
end