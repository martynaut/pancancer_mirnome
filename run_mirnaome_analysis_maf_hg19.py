import click
from maf_files_processing import extract_mirnas_maf, liftover_mutations, dist_occur


@click.command()
@click.argument('input_file')
@click.argument('output_folder')
@click.argument('coordinates_file')
@click.argument('localization_file')
@click.argument('hg38_file')
@click.argument('hg19_file')
@click.option('--from_step', '-s')
def main(input_file,  output_folder, coordinates_file,
         localization_file, hg38_file, hg19_file,
         from_step=''):
    if not from_step:
        from_step = 0
    from_step = int(from_step)
    if from_step <= 1:
        click.echo("Step 1: Extract results for mirnaome")
        extract_mirnas_maf(input_file, output_folder, hg19_file)
        click.echo("miRNAs mutations extracted")
    else:
        click.echo("Skipping step 1")
    if from_step <= 2:
        click.echo("Step 2: Liftover mutations")
        liftover_mutations(input_file, output_folder, hg19_file, hg38_file)
    else:
        click.echo("Skipping step 2")

    if from_step <= 3:
        click.echo("Step 3: Make distinct and occure files")
        dist_occur(output_folder, localization_file, coordinates_file)
    else:
        click.echo("Skipping step 3")

    click.echo("Analysis finished")


if __name__ == "__main__":
    main()
