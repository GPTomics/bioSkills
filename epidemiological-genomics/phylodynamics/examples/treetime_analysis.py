'''Time-scaled phylogenetics with TreeTime'''

from datetime import datetime
import pandas as pd


def convert_to_decimal_year(date_str):
    '''Convert date to decimal year format for TreeTime

    Accepted inputs:
    - ISO date: '2020-06-15'
    - Year-month: '2020-06'
    - Year only: '2020' (treated as midpoint)
    - Already decimal: 2020.5
    '''
    if isinstance(date_str, (int, float)):
        return float(date_str)

    date_str = str(date_str)

    try:
        if '-' in date_str:
            parts = date_str.split('-')
            if len(parts) == 3:
                dt = datetime.strptime(date_str, '%Y-%m-%d')
            elif len(parts) == 2:
                dt = datetime.strptime(date_str + '-15', '%Y-%m-%d')  # Assume mid-month
            else:
                return float(parts[0]) + 0.5
        else:
            return float(date_str)

        year = dt.year
        day_of_year = dt.timetuple().tm_yday
        days_in_year = 366 if (year % 4 == 0) else 365
        return year + (day_of_year - 1) / days_in_year
    except:
        return None


def prepare_dates_file(metadata_df, name_col, date_col, output_path):
    '''Prepare dates file for TreeTime

    Output format:
    name<TAB>date
    sample1<TAB>2020.5
    sample2<TAB>2020.7
    '''
    dates = metadata_df[[name_col, date_col]].copy()
    dates.columns = ['name', 'date']
    dates['date'] = dates['date'].apply(convert_to_decimal_year)
    dates = dates.dropna()

    dates.to_csv(output_path, sep='\t', index=False)
    print(f'Saved {len(dates)} dates to {output_path}')
    return output_path


def run_treetime(tree_file, alignment_file, dates_file, output_dir):
    '''Run TreeTime molecular clock analysis

    Returns dict with key results:
    - root_date: Estimated date of MRCA
    - clock_rate: Substitutions per site per year
    - r_squared: Clock-like behavior (>0.8 = good)
    '''
    import subprocess
    import os

    os.makedirs(output_dir, exist_ok=True)

    cmd = [
        'treetime',
        '--tree', tree_file,
        '--aln', alignment_file,
        '--dates', dates_file,
        '--outdir', output_dir,
        '--coalescent', 'skyline'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f'TreeTime error: {result.stderr}')
        return None

    # Parse results would go here
    return {'output_dir': output_dir}


def interpret_clock_rate(rate, genome_length):
    '''Interpret molecular clock rate

    Typical clock rates by organism:
    - SARS-CoV-2: ~8e-4 subs/site/year (24 subs/year)
    - Influenza A: ~3e-3 subs/site/year
    - HIV: ~2e-3 subs/site/year
    - E. coli: ~1e-7 subs/site/year (5 SNPs/genome/year)
    - M. tuberculosis: ~5e-8 subs/site/year (0.5 SNPs/year)

    For outbreak analysis:
    - Faster rate = more resolution for recent events
    - Slower rate = need longer time span for dating
    '''
    subs_per_year = rate * genome_length

    print('Clock Rate Interpretation')
    print('=' * 40)
    print(f'Rate: {rate:.2e} substitutions/site/year')
    print(f'Genome-wide: ~{subs_per_year:.1f} substitutions/year')

    # Time estimation examples
    for snps in [5, 10, 20, 50]:
        years = snps / (2 * subs_per_year)  # Divide by 2 for pairwise distance
        print(f'  {snps} SNPs between isolates ≈ {years:.1f} years since MRCA')


def estimate_outbreak_origin(root_date, first_sample_date):
    '''Estimate outbreak timing from phylodynamic analysis

    The root date represents the most recent common ancestor (MRCA)
    of the sampled sequences, which may predate the actual outbreak
    origin if early cases were not sampled.
    '''
    time_before_sampling = first_sample_date - root_date

    print(f'\nOutbreak Timing Estimate:')
    print(f'  Root (MRCA) date: {root_date:.2f}')
    print(f'  First sample: {first_sample_date:.2f}')
    print(f'  Time before first sample: {time_before_sampling:.2f} years')
    print(f'  ({time_before_sampling * 365:.0f} days)')


if __name__ == '__main__':
    print('TreeTime Phylodynamics Example')
    print('=' * 50)

    # Example: Interpret a clock rate
    # SARS-CoV-2 example
    rate = 8e-4  # subs/site/year
    genome_length = 29903  # SARS-CoV-2 genome

    interpret_clock_rate(rate, genome_length)

    # Estimate outbreak origin
    root_date = 2019.9  # Late 2019
    first_sample = 2019.98  # Late December 2019

    estimate_outbreak_origin(root_date, first_sample)
