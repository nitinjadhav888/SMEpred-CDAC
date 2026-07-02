# HelixZero-CMS — ICE Cloud Deployment Guide

## For absolute beginners (no Docker knowledge needed)

---

## Part 1: What is Docker? (real-life analogy)

**Imagine this:**

You baked a cake at home. You want to give it to your friend in another city.
- If you just carry the cake in your hands → it will break, melt, or get dirty
- Instead, you put it in a **Tupperware box** → it stays safe and arrives perfect

**Docker is that Tupperware box for software.**

Your SMEPred app needs:
- Python 3.11
- NumPy, LightGBM, FastAPI (libraries)
- The model file (model_b.pkl, 12 MB)
- Configuration files

On your laptop it works because you already installed all these things.
But on a different computer (like ICE Cloud), those things don't exist.

**Docker packs everything into one box (called an "image").**
Then you can send that box to any server and it will work exactly the same.

---

## Part 2: What is ICE Cloud?

ICE Cloud is **C-DAC's own cloud platform** — like Google Cloud or AWS, but made by C-DAC Pune for Indian researchers. It's:
- **Free** (government-subsidized)
- **Has more memory** than Railway (2 GB, 4 GB, or more)
- **Made for exactly this kind of scientific software**

Since you work at C-DAC, you can use ICE Cloud for free.

---

## Part 3: Step-by-step instructions

### Step 0: Check if you have Docker

Open **Windows PowerShell** (press Start, type "PowerShell", open it).

Type:
```powershell
docker --version
```

**If you see something like** `Docker version 24.0.7` → ✅ Docker is installed. Skip to **Step 2**.

**If you see** `'docker' is not recognized` → ❌ Need to install Docker. Go to **Step 1**.

---

### Step 1: Install Docker Desktop (one-time, 5 minutes)

1. Go to: https://www.docker.com/products/docker-desktop/
2. Click **"Download for Windows"**
3. Open the downloaded file (`Docker Desktop Installer.exe`)
4. Click through the installation (keep all default options)
5. When asked, **restart your computer**
6. After restart, Docker will start automatically (you'll see a whale icon in the taskbar)
7. Open PowerShell and type:
   ```powershell
   docker --version
   ```
   You should now see: `Docker version 24.x.x`

✅ Docker is ready.

---

### Step 2: Open the project folder

In PowerShell, go to your SMEPred folder:

```powershell
cd D:\Helixx\smepred
```

---

### Step 3: Build the Docker image (5-10 minutes)

This creates the "Tupperware box" with your app inside.

```powershell
docker build -t helixzero-cms:latest .
```

**What will happen:**
- Docker will download Python 3.11 (takes a minute)
- Then install all the libraries (NumPy, LightGBM, FastAPI, etc.)
- Then copy your code and model files
- Finally, set up the command to start the server

**At the end you should see:**
```
Successfully built abc123...
Successfully tagged helixzero-cms:latest
```

---

### Step 4: Test the container locally (optional but recommended)

Run your container on your own computer first to make sure it works:

```powershell
docker run -d -p 8000:8000 --name helixzero helixzero-cms:latest
```

Now open your browser and go to: http://localhost:8000

You should see the HelixZero-CMS web interface!

If you see it → ✅ Your app works inside Docker.

To stop the test:
```powershell
docker stop helixzero
docker rm helixzero
```

---

### Step 5: Login to ICE Cloud

1. Open your browser and go to: https://icecloud.in/me/dashboard
2. Log in with your C-DAC email and password
3. Look for something called **"Container Registry"** or **"Image Registry"** in the left menu

The registry will give you an address like:
```
icecloud.in/your-username
```

Copy this address — you'll need it.

---

### Step 6: Login to ICE Cloud from PowerShell

In PowerShell, type:

```powershell
docker login icecloud.in
```

It will ask for:
- **Username**: your C-DAC email (e.g., nitin@cdac.in)
- **Password**: your ICE Cloud password

If login succeeds → ✅ You're connected to ICE Cloud.

---

### Step 7: Tag and push your image to ICE Cloud

Now we send your "Tupperware box" to ICE Cloud's storage.

Replace `YOUR_USERNAME` with your actual ICE Cloud username:

```powershell
docker tag helixzero-cms:latest icecloud.in/YOUR_USERNAME/helixzero-cms:latest
docker push icecloud.in/YOUR_USERNAME/helixzero-cms:latest
```

**What this does:**
- `tag` — gives the image a second name that points to ICE Cloud
- `push` — uploads the image to ICE Cloud (like uploading a file to Google Drive)

⏱ This takes 2-3 minutes (the image is about 200 MB).

---

### Step 8: Create a container on ICE Cloud dashboard

1. Go to: https://icecloud.in/me/dashboard
2. Look for **"Containers"** or **"Workloads"** in the left menu → click **"Create Deployment"**
3. Fill in the form:

| Field | What to type |
|-------|-------------|
| **Image name** | `icecloud.in/YOUR_USERNAME/helixzero-cms:latest` |
| **Container name** | `helixzero-cms` |
| **Port** | `8000` |
| **Memory** | `2048` (this is 2 GB — very important!) |
| **CPU** | `1` (or 2 if available) |
| **Environment variable** | `PORT` = `8000` |

4. Click **"Deploy"** or **"Create"**

---

### Step 9: Get your URL

After deployment, ICE Cloud will show you a URL like:

```
https://helixzero-cms-xxxx.icecloud.in
```

Copy this URL and open it in your browser.

You should see the HelixZero-CMS app running! 🎉

---

### Step 10: Test the API

Open a new PowerShell window and test:

```powershell
curl https://helixzero-cms-xxxx.icecloud.in/modifications
```

You should get back a JSON list of all 30 chemical modifications.

---

## Troubleshooting

### "Container keeps crashing / restarting"

**Problem:** Not enough memory.

**Fix:** When creating the container on ICE Cloud, set Memory to at least **2048 MB (2 GB)**. The LightGBM model needs ~150-200 MB in RAM, and during beam search it needs more for the feature matrices.

### "docker: command not found"

**Problem:** Docker is not installed or not running.

**Fix:** 
1. Make sure Docker Desktop is installed (Step 1)
2. Make sure Docker Desktop is running (look for the whale icon in your taskbar)
3. Restart PowerShell and try again

### "Cannot connect to ICE Cloud registry"

**Problem:** Login failed or network issue.

**Fix:**
1. Make sure you can log in at https://icecloud.in/me/dashboard in your browser
2. Try `docker login icecloud.in` again — make sure you type the password correctly
3. If you're on C-DAC office network, you might need to contact the ICE Cloud admin to get registry access

---

## What if I get stuck?

Just ask me. Tell me:
1. Which step you're on
2. What command you typed
3. What error you saw (copy-paste the red text)

I'll help you fix it.

---

## Summary of commands (cheat sheet)

```powershell
# 1. Build the image
cd D:\Helixx\smepred
docker build -t helixzero-cms:latest .

# 2. Test locally
docker run -d -p 8000:8000 --name helixzero helixzero-cms:latest
# Open http://localhost:8000 to test
docker stop helixzero && docker rm helixzero

# 3. Push to ICE Cloud
docker login icecloud.in
docker tag helixzero-cms:latest icecloud.in/YOUR_USERNAME/helixzero-cms:latest
docker push icecloud.in/YOUR_USERNAME/helixzero-cms:latest

# 4. Then go to https://icecloud.in/me/dashboard → Create Container
```
