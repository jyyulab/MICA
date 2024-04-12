FROM ghcr.io/stjudecloud/conda-base:1.1.1 AS builder

RUN conda create -n py376 \
    python==3.7.6 \
    -y \
    && conda clean --all

# Make RUN commands use the new environment:
RUN echo "conda activate py376" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

RUN pip3 install --user --ignore-installed setuptools==57.5.0 \
    && pip3 install --user --ignore-installed mica

FROM debian:10-slim
COPY --from=builder /opt/conda/envs/py376/ /opt/conda/envs/py376/
ENV PATH /opt/conda/envs/py376/bin:$PATH
COPY --from=builder /root/.local/bin/ /root/.local/bin/
COPY --from=builder /root/.local/lib/ /root/.local/lib/
ENV PATH /root/.local/bin:$PATH

ENTRYPOINT [ "mica" ]
