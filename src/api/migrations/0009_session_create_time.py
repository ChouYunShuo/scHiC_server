# Generated by Django 4.1.2 on 2024-08-16 05:52

from django.db import migrations, models
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('api', '0008_session'),
    ]

    operations = [
        migrations.AddField(
            model_name='session',
            name='create_time',
            field=models.DateTimeField(default=django.utils.timezone.now),
        ),
    ]
