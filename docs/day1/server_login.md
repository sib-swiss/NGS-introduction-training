## Material

### Working on the cloud server

!!! warning "Great power comes with great responsibility"
    The cloud server is a temporary instance for this workshop only. Although the computational resources should be more than enough, **it's a basic Ubuntu server, and there are no hard limits on memory or CPU usage.**
    Take therefore into account that great power comes with great responsibility. Overloading it can result in a reboot, cancelling all running calculations.

### Visualisation of results

Many results come in an image (e.g. `.png`, `.jpg`) or `html` format. These can not be viewed directly from the server. Therefore, transfer them to your local PC first (with e.g. [FileZilla](https://filezilla-project.org/) or `scp`), or mount the server with `sshfs` (more info [here](https://www.digitalocean.com/community/tutorials/how-to-use-sshfs-to-mount-remote-file-systems-over-ssh)).

## Exercises

### 1. First login

>:fontawesome-regular-clock: 30 minutes

#### Login to AWS EC2 remote server
You will receive an e-mail shortly before the workshop with a key, username and IP address to login on a cloud server.
Login like this:

```sh
ssh -i path/to/key/key_[USERNAME].pem [USERNAME]@[AWS_IP]
```

#### Setup your favourite editor to work remotely

To directly initiate and modify scripts on the remote server you can use plugins:

* Notepadd++: `NppFTP`
* Atom: `remote-edit-ni`

Setup the connection to the server with the following details:

* protocol: sftp
* username: your username
* hostname: server IP
* port: 22
* authentication/logon type: path to private key file

### 2. A UNIX command line interface (CLI) refresher

>:fontawesome-regular-clock: 30 minutes

Most bioinformatics software are UNIX based and are executed through the CLI. When working with NGS data, it is therefore convenient to improve your knowledge on UNIX. For this course, we need basic understanding of UNIX CLI, so here are some exercises to refresh your memory.

#### Make a new directory

Make a directory called `reads` in your home directory, and make it your current directory.

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
