import smtplib
from email.mime.text import MIMEText

sender = "x" # Replace with your actual email address or use an environment variable for security
receiver = "x" # Replace with the recipient's email address or use an environment variable for security
password = "x" # Replace with your actual password or use an environment variable for security

msg = MIMEText("Your RNA-seq pipeline is DONE!")
msg['Subject'] = "RNA-seq Finished 🎉"
msg['From'] = sender
msg['To'] = receiver

with smtplib.SMTP_SSL('smtp.gmail.com', 465) as server:
    server.login(sender, password)
    server.send_message(msg)

print("Email sent!")