# Installation

To run this software, you need to follow the instructions in [basic requirements](#basic-requirements). This will allow you to run the outward assembly software in the *local* profile (we expect the *local* profile will work for most use-cases).

Briefly, the outward assembly software has two profiles (read our [usage docs](./usage.md#profiles) to learn more about when to use either profile):
* *Local*: Run outward assembly using your local machine.
* *Batch*: Run outward assembly using AWS Batch (via Nextflow) for distributed read searching; other algorithm steps run locally.

In the situation that you want to use the *batch* profile, you will need to follow the [instructions below to install the required software](#optional-batch-profile).

Additionally, if one would like generate high frequency k-mers to use for the outward assembly software, you can follow the [instructions below to install the required software](#optional-kmc).

## Basic requirements

Baseline required dependencies are specified in the conda configuration file `outward_assembly_env.yml`. Use conda or the conda-like package manager of your choice. From the base repository directory: 

```
conda env create -f outward_assembly_env.yml
conda activate OutwardAssembly
# install outward assembly itself in dev profile in the conda env
pip install -e .
```

## (Optional) Batch profile

Using the batch profile requires doing two steps:

1. Install Nextflow and the AWS CLI
2. Configuring AWS Batch

### Nextflow installation

#### 1. Install Nextflow 
To install Nextflow, you can follow the instructions [here](https://www.nextflow.io/docs/latest/getstarted.html)

#### 2. Install AWS CLI
To install the AWS CLI, you can follow the instructions [here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

#### 3. Install Docker
To install Docker, you can follow the instructions [here](https://docs.docker.com/engine/install/).

#### 4. Configure AWS & Docker

[Configure AWS access](https://www.nextflow.io/docs/latest/aws.html) by creating a file at `~/.aws/config` **or** `~/.aws/credentials`, specifying your access key ID and secret access key, e.g.

`~/.aws/config`:
```
[default]
region = us-east-1
output = table
tcp_keepalive = true
aws_access_key_id = <ACCESS_KEY_ID>
aws_secret_access_key = <SECRET_ACCESS_KEY>
```

`~/.aws/credentials`:
```
[default]
aws_access_key_id = <ACCESS_KEY_ID>
aws_secret_access_key = <SECRET_ACCESS_KEY>
```

> [!TIP]
> If you encounter `AccessDenied` errors after doing this, you may also need to export these keys as environment variables before running Nextflow:
>
> ```
> eval "$(aws configure export-credentials --format env)"
> ```

Next, you need to make sure your user is configured to use Docker. To do this, create the `docker` user group and add your current user to it:

```
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

### Setup AWS Batch

To setup AWS Batch, we recommend that you follow the latest [^1] instructions in the [MGS Workflow](https://github.com/naobservatory/mgs-workflow/blob/master/docs/batch.md) [^1], _except_ substitute the following storage configuration in your EC2 launch template:
* Volume type: gp3
* Capacity: 5GB plus whatever space is required for the snapshot used by your launch template
* Throughput: 125MB/s (free tier)
* IOPS: 3000 (free tier)

Though each Batch job runs for just a few seconds, running outward assembly at scale may create hundreds of thousands of Batch jobs. Thus, to control costs, it's important not to provision unnecessarily large EBS volumes. 

[^1]: In the situation that the instructions change, [this link](https://github.com/naobservatory/mgs-workflow/blob/9fe05a5ca9ce7cbc886927788f22c71ff9f26443/docs/batch.md) should be used as reference for the version of the docs used at the time.

### Seqera Tower Credentials
Using the Batch profile requires a Seqera Tower [access token](https://docs.seqera.io/platform-cloud/api/overview). See or create a token in the Seqera [console](https://cloud.seqera.io/tokens).

Occasionally, jobs may fail due to exceeding the [container pull rate limit with Wave](https://docs.seqera.io/wave/api).  To overcome this, you can provide a Tower access token to increase your rate limit by 4x. You can obtain an access token from the Seqera Platform. For more instructions, please refer to this [page](https://github.com/naobservatory/mgs-workflow/blob/master/docs/troubleshooting.md#api-container-errors).

You can provide your user token to Nextflow by setting the token as an environment variable:

```bash
export TOWER_ACCESS_TOKEN=your_tower_access_token
```
When running the pipeline with Nextflow, the token will be automatically used.

Alternatively, the pipeline also supports an optional `tower_token` parameter that can be used to explicitly provide the Tower access token. Passing `tower_token='your_token'` to relevant functions will also configure Nextflow to use it properly.

## (Optional) KMC

The `OutwardAssembly` conda environment does not include the kmer counting tool KMC. Kmer counting with KMC is not required for outward assembly, but you may find it helpful to use KMC to generate a list of high frequency kmers. `kmer_freq_filter.py` contains helper functions for computing high frequency kmers using KMC and creating frequency-filtered subsets of input reads. By default, these functions look for KMC's executable at `outward-assembly/non-conda-deps/kmc3/bin`. KMC binaries are available on [GitHub](https://github.com/refresh-bio/KMC/releases). 
