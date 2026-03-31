# Nonpareil container build notes

This repo uses a Singularity/Apptainer image for Nonpareil at:

- `containers/nonpareil.sif`

## Source image

BioContainer image:

- `quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2`

Note: this repository does not currently publish a `latest` tag, so use an explicit version tag.

## Option A: Pull directly from Docker registry (no definition build)

```bash
cd /storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper
singularity pull containers/nonpareil.sif docker://quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2
```

## Option B: Build from the definition file in this folder

Definition file:

- `containers/definitions/nonpareil.def`

Build command (requires root/sudo or configured fakeroot):

```bash
cd /storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper/containers/definitions
sudo singularity build ../nonpareil.sif nonpareil.def
```

## Quick validation

```bash
cd /storage/biology/projects/miller-lowry/beitner/binning-classification-wrapper
singularity exec containers/nonpareil.sif nonpareil -h
```

## Docker reference (original BioContainers usage)

```bash
docker pull quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2
docker run -it quay.io/biocontainers/nonpareil:3.5.5--r44h077b44d_2 /bin/bash
```
