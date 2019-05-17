import tkFileDialog
import tkFont
from Tkinter import *
import json
import Queue
import threading
from TREEasy import run_TREEasy, dependency_check
import os
import signal
import webbrowser


class TREEasyThread(threading.Thread):
    def __init__(self, fail, finish, **kwargs):
        threading.Thread.__init__(self)
        self.args = kwargs
        self.queue = self.args['message_queue']
        self.fail = fail
        self.finish = finish

    def run(self):
        print("progress start")
        try:
            run_TREEasy(**self.args)
        except Exception, e:
            # self.queue.put("PROCESS FAILED")
            print(e)
            print("progress failed")
            pass


class FileFrame(Frame):
    def __init__(self, master, width, height, title="", isDirectory=False):
        Frame.__init__(self, master=master, width=width, height=height)
        self.pack()
        self.master = master
        self.title = title

        self.label = Label(master=self, text=title, width=20, anchor="e")
        self.label.grid(row=0, column=0)
        self.text = Entry(master=self, width=75)
        self.text.grid(row=0, column=1)
        self.button = Button(master=self, text="select", width=5, command=self.select)
        self.button.grid(row=0, column=2)

        self.isDirectory = isDirectory

    def filename(self):
        return self.text.get()

    def select(self):
        if self.isDirectory:
            path = tkFileDialog.askdirectory(title=self.title)
        else:
            path = tkFileDialog.askopenfilename(title=self.title)
        if path != "":
            self.set(path)

    def set(self, path):
        self.text.delete(0, END)
        self.text.insert(0, path)


class ProgressWindow(Toplevel):
    warning_message = "Warning message here! If you still want to proceed with GUI, click the START button."

    def __init__(self, parent, **kwargs):
        Toplevel.__init__(self, parent)
        self.wm_title("Run TREEasy")

        self.kwargs = kwargs

        # Text Field of shell command
        self.cmd_label = Label(master=self, text="Command Line")
        self.cmd_label.pack()

        self.cmd_field = Text(master=self, height=2)
        self.cmd_field.insert(END, self.get_cmd())
        self.cmd_field.configure(state="disabled")
        self.cmd_field.pack()
        self.grab_set()

        # Progress Text
        self.progress_label = Label(master=self, text="Progress")
        self.progress_label.pack()
        self.progress_field = Text(master=self, height=10)
        self.progress_field.insert(0.0, self.warning_message)
        self.progress_field.pack()
        self.progress_field.configure(state="disabled")

        # Start button
        self.start_button = Button(master=self, text="START", command=self.start)
        self.start_button.pack()
        self.started = False

        # setting up a message queue
        self.message_queue = Queue.Queue()

        # on close
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def get_cmd(self):
        cmd = "python TREEasy.py -d {} -s {} -g {} -b {} -r {} -n {} -k {} -t {} -c {}".format(
            self.kwargs["directory"], self.kwargs["species"], self.kwargs["genes"],
            self.kwargs["bootstrap"], self.kwargs["root"], self.kwargs["network"],
            self.kwargs["cross"], self.kwargs["thread"], self.kwargs["type"]
        )
        print(cmd)
        return cmd

    def start(self):
        self.started = True
        self.start_button.destroy()
        monitor = threading.Thread(target=self.start_monitor)
        monitor.start()
        pass

    def run(self):
        config = self.kwargs
        self.clear()
        self.append("TREEasy started!")
        try:
            run_TREEasy(seq_file=config["directory"], roottaxon=config["root"], type=config["type"],
                        species_namefile=config["species"], gene_namefile=config["genes"],
                        boot_value=config["bootstrap"], Net_num=config["network"],
                        cross_value=config["cross"], thread_number=config["thread"],
                        message_queue=self.message_queue)
            self.message_queue.put("TREEasy finished")
        except Exception, e:
            print(e.message)
            print(Exception.message)
            self.append("TREEasy failed! Please check your parameters.")
            self.message_queue.put("TREEasy Failed")

    def start_monitor(self):
        thread = threading.Thread(target=self.run)
        thread.start()
        while True:
            msg = self.message_queue.get(block=True)
            if msg == "TREEasy Failed" or msg == "TREEasy finished":
                break
            else:
                self.append(msg)
        thread.join()

    def append(self, message):
        self.progress_field.configure(state="normal")
        self.progress_field.insert(END, message + '\n')
        self.progress_field.configure(state="disabled")
        self.progress_field.update()

    def clear(self):
        self.progress_field.configure(state="normal")
        self.progress_field.delete(0.0, END)
        self.progress_field.configure(state="disabled")
        self.progress_field.update()

    def on_close(self):
        print("on close")
        if self.started:
            os.system("killall -9 iqtree")
            os.kill(os.getpid(), signal.SIGKILL)
        self.destroy()
        pass


class TREEasyApp(Tk):
    width = 800

    def __init__(self):
        Tk.__init__(self)
        self.wm_title("TREEasy")
        # default font
        default_font = tkFont.nametofont("TkDefaultFont")
        default_font = tkFont.Font(font=("TimesNewRoman", 12, "bold"))
        default_font.configure(size=10)
        self.option_add("*Font", default_font)
        # menubar
        self.menubar = Menu(self)

        self.filemenu = Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Example - Simulated", command=self.example_simulated)
        self.filemenu.add_command(label="Example - Acropora", command=self.example_acropora)
        self.menubar.add_cascade(label="Examples", menu=self.filemenu)
        self.helpmenu = Menu(self.menubar, tearoff=0)
        self.helpmenu.add_command(label="Help", command=self.help)
        self.menubar.add_cascade(label="Help", menu=self.helpmenu)
        self.config(menu=self.menubar)

        # frame for file / directories
        self.directory_selection = FileFrame(master=self, width=self.width, height=100, title="Gene Locus Directory",
                                             isDirectory=True)

        space = Frame(master=self, height=20).pack()
        self.gene_selection = FileFrame(master=self, width=self.width, height=100, title="Gene Name File")
        space = Frame(master=self, height=20).pack()
        self.species_selection = FileFrame(master=self, width=self.width, height=100, title="Species Name File")
        space = Frame(master=self, height=20).pack()

        # root taxons
        self.parameterFrame = Frame(master=self)
        self.parameterFrame.pack()
        row = 0

        self.root_label = Label(master=self.parameterFrame, width=17, text="Root Taxon(s)", anchor="e")
        self.root_label.grid(row=row, column=0, sticky="E")
        self.root_entry = Entry(master=self.parameterFrame, width=79)
        self.root_entry.grid(row=row, column=1, columnspan=4, sticky="EW")
        row += 1

        space = Frame(master=self.parameterFrame, height=20).grid(row=row)
        row += 1
        self.cross_label = Label(master=self.parameterFrame, width=17, text="Cross Value", anchor="e")
        self.cross_label.grid(row=row, column=0, sticky="E")
        self.cross_entry = Entry(master=self.parameterFrame, width=20)
        self.cross_entry.grid(row=row, column=1, sticky="EW")
        self.net_num = Label(master=self.parameterFrame, width=38, text="Maximum Network Number", anchor="e")
        self.net_num.grid(row=row, column=2, columnspan=2, sticky="E")
        self.net_num_slider = Scale(self.parameterFrame, from_=1, to=5, orient=HORIZONTAL)
        self.net_num_slider.grid(row=row, column=4, sticky="EW")
        row += 1

        space = Frame(master=self.parameterFrame, height=20).grid(row=row)
        row += 1
        self.bootstrap_label = Label(master=self.parameterFrame, width=17, text="Bootstrap Value", anchor="e")
        self.bootstrap_label.grid(row=row, column=0, sticky="E")
        self.bootstrap_scale = Scale(self.parameterFrame, from_=0, to=100, orient=HORIZONTAL)
        self.bootstrap_scale.grid(row=row, column=1, columnspan=4, sticky="EW")
        row += 1

        space = Frame(master=self.parameterFrame, height=20).grid(row=row)
        row += 1
        self.type_label = Label(master=self.parameterFrame, width=17, text="Gene/Locus Type", anchor="e")
        self.type_label.grid(row=row, column=0, sticky="E")
        self.type_value = StringVar()
        choices = ["CDS", "nonCDS"]
        self.type_value.set(choices[0])
        self.type_select = OptionMenu(self.parameterFrame, self.type_value, *choices)
        self.type_select.grid(row=row, column=1, sticky="EW")
        self.thread_label = Label(master=self.parameterFrame, width=17, text="Threads Number", anchor="e")
        self.thread_label.grid(row=row, column=2, sticky="E")
        # try:
        #     import multiprocessing
        #     n_cpu = multiprocessing.cpu_count()
        # except NotImplementedError:
        #     n_cpu = 8
        n_cpu = 16
        self.thread_scale = Scale(self.parameterFrame, from_=1, to=n_cpu, orient=HORIZONTAL)
        self.thread_scale.grid(row=row, column=3, sticky="EW")

        self.start_button = Button(master=self.parameterFrame, text="START", command=self.start)
        self.start_button.grid(row=row, column=4, sticky="EW")

    def example_simulated(self):
        self.load_config("example_simulated.json")

    def example_acropora(self):
        self.load_config("example_acropora.json")

    def start(self):
        run_window = ProgressWindow(self, **self.get_args())

    help_message = "For help on how to use TREEasy, please go to our GitHub website: \n " \
                   "https://github.com/MaoYafei/TREEasy\n" \
                   "Please cite []"

    def help(self):
        help_window = Toplevel(self)
        help_text = Text(help_window, height=4)
        help_text.insert(0.0, self.help_message)
        help_text.pack()
        website_button = Button(help_window, text="Go to website", command=self.open_website)
        website_button.pack()
        help_window.grab_set()

        pass

    def open_website(self):
        webbrowser.open("https://github.com/MaoYafei/TREEasy", new=0, autoraise=True)

    def get_args(self):
        kwargs = {"directory": self.directory_selection.filename(), "genes": self.gene_selection.filename(),
                  "species": self.species_selection.filename(), "bootstrap": self.bootstrap_scale.get(),
                  "cross": int(self.cross_entry.get()), "root": self.root_entry.get(),
                  "network": self.net_num_slider.get(), "thread": self.thread_scale.get(),
                  "type": self.type_value.get()}
        return kwargs

    def load_config(self, filename):
        with open(filename, "rt") as f:
            config = json.load(f)
        self.directory_selection.set(config["directory"])
        self.species_selection.set(config["species"])
        self.gene_selection.set(config["genes"])

        self.root_entry.delete(0, END)
        self.root_entry.insert(0, config["root"])

        self.bootstrap_scale.set(config["bootstrap"])
        self.net_num_slider.set(config["networks"])
        self.cross_entry.delete(0, END)
        self.cross_entry.insert(0, config["cross"])
        self.thread_scale.set(config["thread"])
        self.type_value.set(config["type"])
        # self.type_select.update()


if __name__ == '__main__':
    app = TREEasyApp()
    app.mainloop()
    print("byebye")
    sys.exit(0)
