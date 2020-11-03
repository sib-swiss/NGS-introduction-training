## Material

### Working on the cloud server

In this part we will set up your computer to work on the remote AWS server. You have received an e-mail shortly before the workshop with a key, username and IP address to login on a cloud server. Before you do this part, you should have installed FileZilla, Atom if you're on Linux/Mac OS and MobaXterm if you're on Windows.

!!! warning "Great power comes with great responsibility"
    The cloud server is a temporary instance for this workshop only. Although the computational resources should be more than enough, **it's a basic Ubuntu server, and there are no hard limits on memory or CPU usage.**
    Take therefore into account that great power comes with great responsibility. Overloading it can result in a reboot, cancelling all running calculations.

### Video tutorials

Below you can find video tutorials to set up FileZilla, atom and MobaXterm to edit and/or transfer remote files.

#### Atom

Atom is a versatile text editor for all major operating systems. For this course, it's the recommended script editor for Linux and Mac OS users. With the third-party package `ftp-remote-edit`, you can remotely edit scripts. The video tutorial explains how to set it up.

<iframe src="https://player.vimeo.com/video/473838666" width="640" height="360" frameborder="0" allow="autoplay; fullscreen" allowfullscreen></iframe>

#### FileZilla

Many results come in an image (e.g. `.png`, `.jpg`) or `html` format. These can not be viewed directly from the server. Also, for this course, files loaded in IGV need to be on your local computer. You can easily transfer files between your local PC and the remote host with [FileZilla](https://filezilla-project.org/). The video tutorial explains how to set it up.

<iframe src="https://player.vimeo.com/video/473838726" width="640" height="360" frameborder="0" allow="autoplay; fullscreen" allowfullscreen></iframe>

#### MobaXterm

MobaXterm is an SSH client for Windows. Use this to connect to the remote host and edit remote scripts if you're on Windows. The video tutorial explains how to set it up.

<iframe src="https://player.vimeo.com/video/473838657" width="640" height="360" frameborder="0" allow="autoplay; fullscreen" allowfullscreen></iframe>

## Exercises

### 1. First login

>:fontawesome-regular-clock: 30 minutes

#### Login to AWS EC2 remote server

Use the video tutorials and the information below to log in and set up a remote script editor.

=== "mac OS/Linux"

    Open a terminal and login like this:

    ```sh
    ssh -i path/to/key/key_[USERNAME].pem [USERNAME]@[AWS_IP]
    ```

    !!! warning
        change `path/to/key` to the actual path where you have put the key file.

    #### Setup your favourite editor to work remotely

    To directly initiate and modify scripts on the remote server you can use the Atom plugin `ftp-remote-edit`

    In general, setup the connection to the server with the following details:

    * protocol: sftp
    * username: your username
    * hostname: server IP
    * port: 22
    * authentication/logon type: path to private key file

    Tutorials are found above in the [video tutorial to set up Atom](#Atom).


=== "Windows"
    If you are using MobaXterm on windows, you will automatically login to the remote server once you've started the SSH session. Follow the [video tutorial on MobaXterm](#mobaxterm) to set up an SSH session.

    These are the general settings you should take into account:

    * protocol: sftp
    * username: your username
    * hostname: server IP
    * port: 22
    * authentication/logon type: path to private key file

#### Initiate conda

To make use of the pre-installed software with conda, we need to initiate it first. Login to the server and run:

```sh
/opt/miniconda3/bin/conda init
exec bash
```

To load the environment with the required software packages, run:

```sh
conda activate ngs
```

!!! note "Activating the environment"
    You will need to activate the ngs environment each time you login.

### 2. A UNIX command line interface (CLI) refresher

>:fontawesome-regular-clock: 30 minutes

Most bioinformatics software are UNIX based and are executed through the CLI. When working with NGS data, it is therefore convenient to improve your knowledge on UNIX. For this course, we need basic understanding of UNIX CLI, so here are some exercises to refresh your memory.

#### Make a new directory

Login to the server and use the command line to make a directory called `reads` in your home directory, and make it your current directory.

??? done "Answer"
    ```sh
    cd
    mkdir reads
    cd reads
    ```

#### File permissions

Generate an empty file in your newly made directory like this:

```sh
cd ~/reads
touch new_file.txt
```

Can you use this file as an executable script?

??? done "Answer"
    Use `ls` to get the file permissions:

    ```sh
    ls -l
    ```

    This will give you something like:

    ```
    -rw-r--r--  1 geertvangeest  staff  0 Sep 21 13:38 new_file.txt
    ```

    In any of the user, group or other permissions, there is no `x` for execute. This means it is not executable

#### Redirection: `>` and `|`

In `/home` you'll find a directory for all users with an account. Write all usernames to a file called `usernames.txt`.

??? done "Answer"
    ```sh
    ls /home > ~/usernames.txt
    ```

The command `wc -l` counts the number of lines, and can read from stdin. Make a one-liner with a pipe `|` symbol to find out how many users have an account on this server.

??? done "Answer"
    ```sh
    ls /home | wc -l
    ```

#### Variables

Store `usernames.txt` as variable (like this: `VAR=variable`), and use `wc -l` on that variable to count the number of lines in the file.

??? done "Answer"
    ```sh
    FILE=usernames.txt
    wc -l $FILE
    ```

#### shell scripts

Make a shell script that automatically counts the current number of account-holders.

??? done "Answer"
    Make a script called e.g. `current_accounts.sh`:
    ```sh
    #!/usr/bin/env bash
    cd /home
    ls | wc -l
    ```

### 3. Detaching a job

>:fontawesome-regular-clock: 20 minutes

On this server, there is no job scheduler, so everything is run directly from the command line. That means that if a process is running, the command line will be busy, and the job will be killed upon logout. To circumvent this, there are several methods to 'detach' the screen or prevent a 'hangup signal' of a job runnig in the background that will terminate your running job.
The software `screen` or `tmux` can be used to detach your screen, and all messages to stderr or stdout (if not redirected) will be printed to the (detached) console. Use those if you're comfortable with them.

Another, more basic, program to prevent the 'hangup signal' is `nohup`. Use it like so:

```sh
nohup [YOUR COMMAND] &
```

!!! note
    Don't forget the `&` after the command. This symbol let's the process run in the background.

So, for running e.g. a shell script this would be:

```sh
nohup script.sh &
```

Anything written to stdout or stderr will be written to the file `nohup.out` in your current working directory.

!!! warning "Don't change your script while running"
    `nohup` runs through your script line-by-line and reads it from disk. If you change your script while running it, it will run the new lines.

Generate a script that waits for 120 seconds (use `sleep`) and prints `I'm done!` to stdout if it's done.

??? done "Answer"
    Script called `wait.sh`:
    ```sh
    #!/usr/bin/env bash
    sleep 120
    echo "I'm done!"
    ```

Run this script in the background without a hangup signal (so, by using `nohup`), logout, and login again, and see if it's still running with `ps x`.

??? done "Answer"
    ```sh
    nohup wait.sh &
    ```
    Remember the PID, logout and login (within 120 seconds.. )

    ```sh
    ps x
    ```
