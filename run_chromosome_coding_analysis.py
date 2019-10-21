import click
from prepare_vcf_files import make_unique_files
from coding_chromosomes_analysis import all_files_processing_coding, merge_data


@click.command()
@click.argument('input_folder')
@click.argument('output_folder')
@click.argument('coordinates_file')
@click.option('--from_step', '-s')
def main(input_folder,  output_folder, coordinates_file, from_step=''):
    if not from_step:
        from_step = 0
    from_step = int(from_step)
    if from_step <= 1:
        click.echo("Step 1: Analysis started")
        make_unique_files(input_folder=input_folder, output_folder=output_folder)
        click.echo("Not unique files combined")
    else:
        click.echo("Skipping step 1")
    if from_step <= 2:
        click.echo("Step 2: Extract results for coding")
        all_files_processing_coding(input_folder, output_folder, coordinates_file)
    else:
        click.echo("Skipping step 2")
    if from_step <= 3:
        click.echo("Step 3: Merging data")
        merge_data(output_folder, coordinates_file)
    else:
        click.echo("No steps performed")
    click.echo("Analysis finished")


if __name__ == "__main__":
    main()
