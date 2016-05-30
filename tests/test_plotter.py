from moca.plotter import create_plot
import os
import pytest

@pytest.mark.mpl_image_compare(baseline_dir='data/images',
                               filename='ENCSR000AKB_PhyloP_1.png')

def join_path(head, leaf):
    return os.path.join(head, leaf)

def test_image():
    base_path = 'tests/data/ENCSR000AKB/'
    meme_file = base_path+'moca_output/meme_out/meme.txt'
    plot_title = 'ENCSR000AKB.sorted'
    oc = 'tests/data/generate_out'
    motif_number = 1
    flank_motif = 5
    sample_score_files = [base_path+'moca_output/fimo_out_1/phylop.mean.txt']
    control_score_files = [base_path+'moca_output/fimo_random_1/phylop.mean.txt']
    plot_titles = ['PhyloP']
    centrimo_dir = base_path + 'moca_output/centrimo_out'
    figs = create_plot(meme_file,
                       plot_title,
                       output_dir=oc,
                       centrimo_dir=centrimo_dir,
                       motif_number=motif_number,
                       flank_length=flank_motif,
                       sample_score_files=sample_score_files,
                       control_score_files=control_score_files,
                       reg_plot_titles=plot_titles,
                       annotate=None,
                       save=False)
    return figs[0]


