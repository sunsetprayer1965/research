# Base image
FROM python:3.8-slim


# Install system dependencies required for scientific libraries
RUN apt-get update && apt-get install -y \
    build-essential \
    libatlas-base-dev \
    liblapack-dev \
    gfortran \
    gcc \
    make \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Install required Python libraries using pip
RUN pip install --upgrade pip
RUN pip install numpy==1.14.5
RUN pip install scipy scikit-learn pandas matplotlib iced pastis openpyxl

# Copy the local requirements.txt (if needed)
COPY . .

# Install additional dependencies from requirements.txt
RUN pip install -r requirements.txt
RUN make

# Default command
CMD ["python"]
## docker build -t pastis .
## docker r
## docker run -it --rm -v ${PWD}:/app pastis python /app/structure.py
## docker run -it --rm -v ${PWD}:/app pastis pastis-mds . 
## Get-Content ".\runs\EcoRI+Msel\logs\pastis_pm2.log" -Wait