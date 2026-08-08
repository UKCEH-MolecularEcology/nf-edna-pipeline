process HONEYPI_DOWNLOAD_DB {
    tag "download_${its_region}_db"
    label 'process_low'

    container 'python:3.10-slim'

    storeDir "${params.honeypi.rdp_db_cache}"

    input:
    val db_url
    val its_region

    output:
    path "rdp_db_${its_region}/", emit: db_dir

    script:
    """
    #!/usr/bin/env python3
    import urllib.request
    import tarfile
    import os
    import sys

    url     = "${db_url}"
    region  = "${its_region}"
    out_dir = f"rdp_db_{region}"
    archive = f"rdp_db_{region}.tar.gz"

    print(f"Downloading RDP training data for {region} from:", url, flush=True)

    def _progress(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0 and downloaded % (50 * 1024 * 1024) < block_size:
            pct = min(100, downloaded / total_size * 100)
            print(f"  {pct:.0f}% ({downloaded // 1024 // 1024} MB)", flush=True)

    try:
        urllib.request.urlretrieve(url, archive, reporthook=_progress)
    except Exception as e:
        sys.exit(f"Download failed: {e}\\nURL: {url}\\nProvide a local database via --honeypi.rdp_db_dir instead.")

    print("Extracting archive...", flush=True)
    with tarfile.open(archive) as tf:
        tf.extractall(".")

    extracted = [d for d in os.listdir(".") if os.path.isdir(d) and d != out_dir]
    if not extracted:
        sys.exit("Extraction produced no directories.")
    if len(extracted) > 1:
        extracted.sort(key=lambda d: (
            any(k in d.lower() for k in ('gweon', 'rdp', 'train')), len(d)
        ), reverse=True)
    os.rename(extracted[0], out_dir)
    os.remove(archive)

    props = [f for f in os.listdir(out_dir) if f.endswith('.properties')]
    if not props:
        sys.exit(f"No .properties file found in {out_dir}. Check archive contents.")
    print(f"Database ready: {out_dir}/{props[0]}", flush=True)

    with open(os.path.join(out_dir, "versions.yml"), "w") as fh:
        fh.write('"HONEYPI_DOWNLOAD_DB":\\n')
        fh.write(f'    database_url: {url}\\n')
    """

    stub:
    """
    mkdir -p "rdp_db_${its_region}"
    touch "rdp_db_${its_region}/Gweon-${its_region}.properties"
    """
}
