FROM python:3.6.4

MAINTAINER "Matt Lyon" matt.lyon@bristol.ac.uk

# copy flask app to container
COPY ./requirements.txt /app/
WORKDIR /app

# install python dependencies
RUN pip install --upgrade pip
RUN pip install -r requirements.txt

RUN wget -O /bin/watcher.py https://raw.githubusercontent.com/MRCIEU/bgc-upload-orchestrator/master/watcher.py?token=AB1fTMa-e8fsZrhTfgrH2VEnYtpvjtCBks5cgZIdwA%3D%3D && chmod 775 /bin/watcher.py

# Path
ENV PATH /app:$PATH

CMD tail -f /dev/null
# launch app
# CMD ["python", "main.py"]
