# To deploy these completions you must
# make an alias for relate by adding to ~/.zshrc:
#
#     alias relate='/path/to/relate/bin/Relate'
#
# and then copy and paste this script to the bottom of ~/.zshrc with:
#
#     cat .relate_completion.zsh >> ~/.zshrc

autoload -U compinit; compinit
function _relate(){
    _describe 'command' "(\
        '--mode:Choose which part of the algorithm to run.' \
        '--haps:Filename of the haps file (Output file format of Shapeit).'\
        '--sample:Filename of the sample file (Output file format of Shapeit).'\
        '--map:Genetic map.'\
        '--mutation_rate:Mutation rate.'\
        '--effectiveN:Effective population size.'\
        '--output:Filename of outputu without file extension.'\
        '--dist:Optional but recommended. Distance in BP \nbetween SNPs. Can be generated using RelateFileFormats. If unspecified, distances in haps are used.'\
        '--annot:Optional. Filename of file containing additional annotation of snps. Can be generated using RelateFileFormats.'\
        '--memory:Optional. Approximate memory allowance in GB for storing distance matrices. Default is 5GB.'\
        '--sample_ages:Optional. Filename of file containing sample ages (one per line).'\
        '--chunk_index:Optional. Index of chunk. (Use when running parts of the algorithm on an individual chunk).'\
        '--first_section:Optional. Index of first section to infer. (Use when running parts of algorithm on an individual chunk.'\
        '--last_section:Optional. Index of last section to infer. (Use when running parts of algorithm on an individual chunk.'\
        '--coal:Optional. Filename of file containing coalescent rates. If specified, it will overwrite --effectiveN.'\
        '--transversion:Only use transversion for bl estimation.'\
        '--input:Filename of input'\
        '--painting:Optional. Copying and transition parameters in chromosome painting algorithm. Format: theta, rho. Default: 0.025, 1.'\
        '--seed:Optional. Seed for MCMC im branch lengths estimation.'\
    )"
        
}

compdef _relate relate
setopt complete_aliases
