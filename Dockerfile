FROM python:3.6.4

MAINTAINER "Matt Lyon" matt.lyon@bristol.ac.uk

# copy flask app to container
COPY ./requirements.txt /app/
WORKDIR /app

# install python dependencies
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

# launch app
CMD ["python", "main.py"]
