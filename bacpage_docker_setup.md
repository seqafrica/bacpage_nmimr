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

# Updating the BacPage Docker image 
1. Pull the latest changes from the repository
```bash
git pull
```
2. Rebuild the Docker image
```bash 
docker build -t bacpage:latest .
```

# Updating BacPage
1. Navigate to the directory where you cloned the BacPage repository on the command line
```bash
cd bacpage/
```
