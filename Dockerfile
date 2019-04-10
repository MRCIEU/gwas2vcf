FROM python:3.6.4

MAINTAINER "Matt Lyon" matt.lyon@bristol.ac.uk

# copy flask app to container
COPY ./requirements.txt /app/
WORKDIR /app

# install python dependencies
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

ADD watcher.py /home/bin
RUN chmod 775 /home/bin/watcher.py

# Path
ENV PATH /app:$PATH

CMD tail -f /dev/null
# launch app
# CMD ["python", "main.py"]
