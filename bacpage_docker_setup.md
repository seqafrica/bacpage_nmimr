# Running BacPage with Docker

If you prefer to run BacPage in a Docker container, follow these steps:

#### Building the Docker Image

1. Clone the BacPage repository:
   ```bash
   git clone https://github.com/seqafrica/bacpage_nmimr.git
   cd bacpage_nmimr
   ```

## Build the Docker image
```bash
docker build -t bacpage:latest .
```
## Running BacPage with Docker
- To run BacPage, use the following command to start a Docker container:
```bash
docker run --rm -it -v /path/to/your/data:/data bacpage:latest bacpage -h
```
