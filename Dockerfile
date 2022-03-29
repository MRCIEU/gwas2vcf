FROM python:3.8

MAINTAINER "Matt Lyon" matt.lyon@bristol.ac.uk

# copy flask app to container
COPY . /app
WORKDIR /app

# install python dependencies
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN pip install git+https://github.com/bioinformed/vgraph@v1.4.0#egg=vgraph

# launch app
CMD ["python", "main.py"]
