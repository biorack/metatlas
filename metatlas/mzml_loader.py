import pymzml
import os
import pwd
import datetime
import tables

DEBUG = False


class Spectrum(tables.IsDescription):
    mz = tables.Float32Col()
    scan_time = tables.Float32Col()
    i = tables.Float32Col()
    polarity = tables.UInt8Col()
    ms_level = tables.UInt8Col()
    precursor_MZ = tables.Float32Col()
    precursor_intensity = tables.Float32Col()
    collision_energy = tables.Float32Col()


def mzml_to_hdf(in_file_name, out_file_name=None):
    """Converts in_file (mzml) to binary and stores it in out_file
    """
    if not out_file_name:
        out_file_name = in_file_name.replace('.mzML', '.h5')

    out_file = tables.open_file(out_file_name, "w")
    table = out_file.create_table('/', 'spectra', description=Spectrum)
    if DEBUG:
        print("STATUS: Converting %s to %s (mzML to HDF)" %
              (in_file_name, out_file_name))

    # Extra accessions for pymzml to read
    extraAccessions = [
        ('MS:1000016', ['value', 'unitName']),  # scan start time
        ('MS:1000129', ['name']),  # negative scan
        ('MS:1000130', ['name']),  # positive scan
        ('MS:1000744', ['name', 'value']),  # selected ion m/z
        ('MS:1000042', ['name', 'value']),  # peak intensity
        ('MS:1000045', ['name', 'value']),  # collision energy
    ]

    mzml_reader = pymzml.run.Reader(in_file_name,
                                    extraAccessions=extraAccessions)

    iteration = 0

    for spectrum in mzml_reader:
        try:
            polarity = 1 if 'positive scan' in spectrum.keys() else 0
            ms_level = spectrum['ms level']
            # check if scan start time exists. some thermo spectra
            #  are missing this value
            if 'scan start time' in spectrum.keys():
                scan_time = spectrum['scan start time'][0]
            else:
                scan_time = spectrum['MS:1000016'][0]
        except KeyError:
            continue

        for mz, i in spectrum.peaks:
            precursor_MZ = 0.0
            precursor_intensity = 0.0
            collision_energy = 0.0

            if ms_level == 2:
                collision_energy = spectrum['collision energy'][1]
                if 'peak intensity' in spectrum.keys():
                    precursor_intensity = spectrum['peak intensity'][1]
                else:
                    precursor_intensity = 0.0
                precursor_MZ = spectrum['selected ion m/z'][0]

            table.append([((int(mz), float(scan_time), float(i),
                            int(polarity), int(ms_level),
                            float(precursor_MZ), float(precursor_intensity),
                            float(collision_energy)))])
            if not (iteration % 1024):
                table.flush()

    if DEBUG:
        print('Saving file')
    out_file.set_node_attr('/', "upload_date", datetime.datetime.utcnow())
    out_file.set_node_attr('/', "uploaded_by",
                           pwd.getpwuid(os.getuid())[0])
    out_file.close()
    os.chmod(out_file_name, 770)
    if DEBUG:
        print("STATUS: Finished mzML to HDF conversion")


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Load mzml files to HDF")
    parser.add_argument("-o", "--output", nargs=1, type=str,
                        help="Output file name", required=False)
    parser.add_argument("-i", "--input", type=str, nargs=1,
                        help="Input mzML file", required=True)
    parser.add_argument("-d", "--debug", help="Sets debug mode",
                        action="store_true")

    args = parser.parse_args()
    kwargs = {}
    if args.output:
        kwargs['output_name'] = args.output[0]
    if args.input:
        kwargs['input_name'] = args.array[0]

    # Toggles debug mode base on --debug flag
    DEBUG = args.debug

    mzml_to_hdf(args.input, args.output)
